#!/usr/bin/env python3
# =============================================================================
# PUNTS DE CALOR PER SATÈL·LIT — NASA FIRMS — La TdG
# Portat de FIRMSLaTdG.gs (projecte Correu d'Apps Script) a GitHub el 21.10.2026,
# per estalviar quota de Google. Fa el mateix, amb els mateixos textos:
#   · consulta FIRMS (VIIRS S-NPP, NOAA-20, NOAA-21 i MODIS) a 30 km de la Bisbal;
#   · avisa el bot d'infos de les deteccions noves (les de confiança baixa i els
#     punts ignorats, no);
#   · marca els punts recurrents (possible font fixa);
#   · per als focus seriosos, encarrega a satellit.yml (aquest mateix repo) la
#     zona cremada amb Sentinel-2.
# La memòria (deteccions vistes, històric de punts i encàrrecs) es desa a la
# branca «firms-estat». La primera execució només marca les deteccions que ja
# hi ha, sense enviar res.
#
# Variables d'entorn: FIRMS_MAP_KEY, TELEGRAM_TOKEN_INFOS, TELEGRAM_CHAT_INFOS,
# GH_TOKEN, REPO.  Ús: python3 firms.py <carpeta de l'estat>
# =============================================================================

import csv
import datetime as dt
import io
import json
import math
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from zoneinfo import ZoneInfo

TZ = ZoneInfo("Europe/Madrid")
URL_API = "https://firms.modaps.eosdis.nasa.gov/api/area/csv/"
URL_MAPA = "https://firms.modaps.eosdis.nasa.gov/map/#d:24hrs;@3.04,41.96,10z"
CENTRE = {"nom": "la Bisbal", "lat": 41.9597, "lon": 3.0386}
RADI_KM = 30
DIES_CONSULTA = 2
FONTS = [("VIIRS_SNPP_NRT", "VIIRS (S-NPP)"), ("VIIRS_NOAA20_NRT", "VIIRS (NOAA-20)"),
         ("VIIRS_NOAA21_NRT", "VIIRS (NOAA-21)"), ("MODIS_NRT", "MODIS")]
MODIS_CONFIANCA_MINIMA = 30
DIES_RECURRENT = 3
# Punts que no s'avisen mai (fonts industrials conegudes). Exemple:
#   {"nom": "Ceràmica X", "lat": 41.9601, "lon": 3.0412, "radiKm": 0.5}
PUNTS_IGNORATS = []
SATELLIT_MIN_DETECCIONS = 2
SATELLIT_MIN_MW = 8


def peticio(url, metode="GET", dades=None, capcaleres=None):
    req = urllib.request.Request(url, data=dades, headers={"User-Agent": "LaTdG-firms/1.0", **(capcaleres or {})},
                                 method=metode)
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()
    except Exception as e:
        return 0, str(e).encode()


def dist_km(lat1, lon1, lat2, lon2):
    rad = math.pi / 180
    a = (math.sin((lat2 - lat1) * rad / 2) ** 2 +
         math.cos(lat1 * rad) * math.cos(lat2 * rad) * math.sin((lon2 - lon1) * rad / 2) ** 2)
    return 6371 * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def rumb(lat1, lon1, lat2, lon2):
    rad = math.pi / 180
    y = math.sin((lon2 - lon1) * rad) * math.cos(lat2 * rad)
    x = math.cos(lat1 * rad) * math.sin(lat2 * rad) - math.sin(lat1 * rad) * math.cos(lat2 * rad) * math.cos((lon2 - lon1) * rad)
    graus = (math.atan2(y, x) / rad + 360) % 360
    return ["al nord", "al nord-est", "a l'est", "al sud-est", "al sud", "al sud-oest", "a l'oest",
            "al nord-oest"][int(math.floor(graus / 45 + 0.5)) % 8]


def cella(lat, lon):
    return f"{round(lat / 0.005) * 0.005:.3f},{round(lon / 0.006) * 0.006:.3f}"


def ignorat(d):
    return next((p for p in PUNTS_IGNORATS if dist_km(p["lat"], p["lon"], d["lat"], d["lon"]) <= p.get("radiKm", 0.5)), None)


def rectangle():
    dlat = RADI_KM / 111.32
    dlon = RADI_KM / (111.32 * math.cos(CENTRE["lat"] * math.pi / 180))
    return ",".join(f"{n:.4f}" for n in (CENTRE["lon"] - dlon, CENTRE["lat"] - dlat, CENTRE["lon"] + dlon, CENTRE["lat"] + dlat))


def consultar_font(clau, codi, nom, dies):
    codi_http, cos = peticio(f"{URL_API}{clau}/{codi}/{rectangle()}/{dies}")
    text = cos.decode("utf-8", errors="replace").strip()
    if codi_http != 200:
        raise RuntimeError(f"HTTP {codi_http} — {text[:200]}")
    if not text.lower().startswith("latitude"):
        raise RuntimeError(text[:200] or "resposta buida")
    sortida = []
    for f in csv.DictReader(io.StringIO(text)):
        f = {k.strip().lower(): (v or "").strip() for k, v in f.items()}
        try:
            lat, lon = float(f["latitude"]), float(f["longitude"])
            hhmm = f.get("acq_time", "").zfill(4)
            quan = dt.datetime.fromisoformat(f"{f['acq_date']}T{hhmm[:2]}:{hhmm[2:]}:00+00:00")
        except (KeyError, ValueError):
            continue
        conf = f.get("confidence", "").lower()
        if conf.isdigit():
            confianca, acceptada = conf + " %", int(conf) >= MODIS_CONFIANCA_MINIMA
        else:
            confianca, acceptada = {"l": "baixa", "n": "normal", "h": "alta"}.get(conf, conf), conf != "l"
        try:
            frp = float(f.get("frp", ""))
        except ValueError:
            frp = None
        km = dist_km(CENTRE["lat"], CENTRE["lon"], lat, lon)
        if km > RADI_KM:
            continue
        sortida.append({
            "id": f"{codi}|{f['acq_date']}|{hhmm}|{lat:.4f}|{lon:.4f}", "font": nom, "lat": lat, "lon": lon,
            "quan": quan.timestamp() * 1000, "dia": f["acq_date"], "confianca": confianca, "acceptada": acceptada,
            "frp": frp, "nit": f.get("daynight", "").upper() == "N", "km": km,
            "rumb": rumb(CENTRE["lat"], CENTRE["lon"], lat, lon),
        })
    return sortida


def consultar(clau, dies):
    totes, correctes = [], 0
    for codi, nom in FONTS:
        try:
            totes += consultar_font(clau, codi, nom, dies)
            correctes += 1
        except Exception as e:
            print(f"⚠️ FIRMS ({nom}): no s'ha pogut consultar ({e}).")
    return sorted(totes, key=lambda d: d["quan"]) if correctes else None


def apuntar_historic(historic, deteccions):
    for d in deteccions:
        c = cella(d["lat"], d["lon"])
        historic.setdefault(c, [])
        if d["dia"] not in historic[c]:
            historic[c].append(d["dia"])
    limit = (dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=30)).strftime("%Y-%m-%d")
    for c in list(historic):
        historic[c] = [x for x in historic[c] if x >= limit]
        if not historic[c]:
            del historic[c]


def dies_punt(historic, d):
    return len(historic.get(cella(d["lat"], d["lon"]), []))


def linia(d, historic):
    hora = dt.datetime.fromtimestamp(d["quan"] / 1000, TZ).strftime("%d.%m.%Y %H.%M") + " h"
    lloc = (f"a menys d'1 km de {CENTRE['nom']}" if d["km"] < 1
            else f"a {int(math.floor(d['km'] + 0.5))} km {d['rumb']} de {CENTRE['nom']}")
    dies = dies_punt(historic, d)
    return ("• " + hora + (" (de nit)" if d["nit"] else "") + " · " + lloc +
            f"\n  {d['font']} · confiança {d['confianca']}" +
            (f" · {d['frp']:.1f}".replace(".", ",") + " MW" if d["frp"] is not None else "") +
            (f"\n  🔁 Punt recurrent ({dies} dies en un mes): possible font fixa" if dies >= DIES_RECURRENT else "") +
            f"\n  📍 https://www.google.com/maps?q={d['lat']:.5f},{d['lon']:.5f}")


def missatges(deteccions, historic):
    n = len(deteccions)
    cap = (f"🔥 PUNTS DE CALOR DETECTATS PER SATÈL·LIT (NASA FIRMS)\n{n} " + ("detecció" if n == 1 else "deteccions") +
           f" a menys de {RADI_KM} km de {CENTRE['nom']}")
    peu = ("\nℹ️ Dades de satèl·lit amb unes 3 h de retard. Pot ser un incendi, una crema o una font industrial: "
           "cal comprovar-ho.\n🗺️ " + URL_MAPA)
    sortida, actual = [], cap
    for d in deteccions:
        l = "\n" + linia(d, historic)
        if len(actual) + len(l) + len(peu) > 3800:
            sortida.append(actual)
            actual = cap + " (continuació)"
        actual += l
    sortida.append(actual + "\n" + peu)
    return sortida


def enviar(text):
    token, xat = os.environ.get("TELEGRAM_TOKEN_INFOS", ""), os.environ.get("TELEGRAM_CHAT_INFOS", "")
    if not token or not xat:
        print("🧪 [sense secrets de Telegram, no s'envia]\n" + text)
        return False
    dades = urllib.parse.urlencode({"chat_id": xat, "text": text, "disable_web_page_preview": "false"}).encode()
    codi, cos = peticio(f"https://api.telegram.org/bot{token}/sendMessage", "POST", dades,
                        {"Content-Type": "application/x-www-form-urlencoded"})
    try:
        if codi == 200 and json.loads(cos).get("ok"):
            return True
    except Exception:
        pass
    print(f"❌ Telegram: HTTP {codi} {cos[:200]!r}")
    return False


def encarregar_satellit(deteccions, historic, encarregats):
    focus = []
    for d in deteccions:
        f = next((g for g in focus if dist_km(g["lat"], g["lon"], d["lat"], d["lon"]) <= 2), None)
        if f:
            f["punts"].append(d)
            f["lat"] = sum(p["lat"] for p in f["punts"]) / len(f["punts"])
            f["lon"] = sum(p["lon"] for p in f["punts"]) / len(f["punts"])
        else:
            focus.append({"lat": d["lat"], "lon": d["lon"], "punts": [d]})

    for f in focus:
        max_mw = max(p["frp"] or 0 for p in f["punts"])
        recurrent = all(dies_punt(historic, p) >= DIES_RECURRENT for p in f["punts"])
        if not ((len(f["punts"]) >= SATELLIT_MIN_DETECCIONS or max_mw >= SATELLIT_MIN_MW) and not recurrent):
            continue
        dia = sorted(p["dia"] for p in f["punts"])[0]
        clau = cella(f["lat"], f["lon"]) + "|" + dia
        if clau in encarregats:
            continue
        km = dist_km(CENTRE["lat"], CENTRE["lon"], f["lat"], f["lon"])
        nom = (f"Focus a {int(math.floor(km + 0.5))} km {rumb(CENTRE['lat'], CENTRE['lon'], f['lat'], f['lon'])} "
               f"de {CENTRE['nom']} ({dia[8:10]}.{dia[5:7]}.{dia[:4]})")
        cos = json.dumps({"ref": "main", "inputs": {
            "mode": "incendi", "lat": f"{f['lat']:.5f}", "lon": f"{f['lon']:.5f}", "data": dia, "nom": nom, "radi": "3",
            "punts": json.dumps([[round(p["lat"], 5), round(p["lon"], 5)] for p in f["punts"][:20]])}}).encode()
        codi, resp = peticio(f"https://api.github.com/repos/{os.environ.get('REPO', 'latdg/satellit')}/actions/workflows/"
                             "satellit.yml/dispatches", "POST", cos,
                             {"Authorization": "Bearer " + os.environ.get("GH_TOKEN", ""),
                              "Accept": "application/vnd.github+json", "Content-Type": "application/json"})
        if codi == 204:
            encarregats[clau] = time.time() * 1000
            print("🛰️ FIRMS: encarregada l'anàlisi de Sentinel-2 — " + nom)
            enviar(f"🛰️ He encarregat l'anàlisi de la zona cremada amb Sentinel-2 ({nom}). "
                   "Arribarà quan hi hagi imatges sense núvols d'abans i de després.")
        else:
            print(f"⚠️ FIRMS: GitHub no ha acceptat l'encàrrec de Sentinel-2 (HTTP {codi}): {resp[:200]!r}")
    limit = time.time() * 1000 - 15 * 86400000
    for k in [k for k, v in encarregats.items() if v < limit]:
        del encarregats[k]


def main(carpeta):
    carpeta.mkdir(parents=True, exist_ok=True)
    fitxer = carpeta / "estat.json"
    try:
        estat = json.loads(fitxer.read_text(encoding="utf-8"))
        primera = False
    except Exception:
        estat, primera = {}, True
    vistes = estat.setdefault("vistes", {})
    historic = estat.setdefault("historic", {})
    encarregats = estat.setdefault("encarregats", {})

    clau = os.environ.get("FIRMS_MAP_KEY", "").strip()
    if not clau:
        raise SystemExit("Falta el secret FIRMS_MAP_KEY.")
    deteccions = consultar(clau, DIES_CONSULTA)
    if deteccions is None:
        raise SystemExit("❌ FIRMS: no s'ha pogut consultar cap satèl·lit.")

    apuntar_historic(historic, deteccions)
    noves = [d for d in deteccions if d["id"] not in vistes]
    avisar = []
    for d in noves:
        if not d["acceptada"]:
            continue
        ign = ignorat(d)
        if ign:
            print(f"FIRMS: detecció a {ign['nom']} (punt ignorat).")
            continue
        avisar.append(d)
    print(f"FIRMS: {len(deteccions)} detecció(ns) a menys de {RADI_KM} km en {DIES_CONSULTA} dies · "
          f"{len(noves)} de noves · {len(avisar)} per avisar.")

    enviat = True
    if primera:
        print("FIRMS: primera revisió; les deteccions actuals es marquen com a vistes sense enviar res.")
    elif avisar:
        enviat = all(enviar(t) for t in missatges(avisar, historic))
        print("✅ FIRMS: avís enviat al bot d'infos." if enviat else "⚠️ FIRMS: Telegram ha fallat; es tornarà a provar.")
        if enviat:
            try:
                encarregar_satellit(avisar, historic, encarregats)
            except Exception as e:
                print(f"⚠️ FIRMS: no s'ha pogut encarregar Sentinel-2 ({e}).")

    if enviat:
        for d in noves:
            vistes[d["id"]] = d["quan"]
        limit = time.time() * 1000 - 4 * 86400000
        for k in [k for k, v in vistes.items() if v < limit]:
            del vistes[k]
    fitxer.write_text(json.dumps(estat, separators=(",", ":")), encoding="utf-8")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "firms-estat"))
