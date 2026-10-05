# =============================================================================
# SATÈL·LIT — La TdG · 05.10.2026
# Repositori públic latdg/satellit (workflow satellit.yml).
#
# MÒDUL 1: ZONES CREMADES AMB SENTINEL-2
#   Per a un incendi (punt i data), busca la imatge sense núvols més propera
#   ABANS de l'incendi i la més recent DESPRÉS, calcula l'índex de crema (dNBR),
#   la superfície cremada i la gravetat, i envia al bot d'infos les dues imatges
#   (la de després, amb la zona cremada pintada) i un resum.
#   Amb data de fi (incendis passats, canvi 05.10.2026), la imatge de després és
#   la PRIMERA sense núvols després que s'apagués (la vegetació encara no ha
#   rebrotat); sense data de fi, la més recent.
#   Si encara no hi ha cap imatge neta de després, l'incendi queda a
#   pendents.json i es torna a provar cada dia (fins a DIES_MAX_PENDENT dies).
#
# Ordres:
#   python3 satellit.py incendi    (dades a les variables LAT, LON, DATA, DATA_FI, NOM, RADI, PUNTS)
#   python3 satellit.py pendents   (revisió diària dels incendis pendents)
#
# Font: Copernicus Data Space (Sentinel Hub), amb CDSE_CLIENT_ID i CDSE_CLIENT_SECRET.
# =============================================================================
import io, json, math, os, sys, time
from datetime import date, datetime, timedelta, timezone

import numpy as np
import requests

CLIENT_ID = os.environ.get("CDSE_CLIENT_ID", "").strip()
CLIENT_SECRET = os.environ.get("CDSE_CLIENT_SECRET", "").strip()
TG_TOKEN = os.environ.get("TELEGRAM_TOKEN_INFOS", "").strip()
TG_CHAT = os.environ.get("TELEGRAM_CHAT_INFOS", "").strip() or "8739529503"

URL_TOKEN = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
URL_CATALEG = "https://sh.dataspace.copernicus.eu/api/v1/catalog/1.0.0/search"
URL_PROCES = "https://sh.dataspace.copernicus.eu/api/v1/process"
CRS_UTM = "http://www.opengis.net/def/crs/EPSG/0/32631"

PENDENTS = "pendents.json"
DIES_ABANS = 45              # finestra per buscar la imatge d'abans
DIES_MAX_PENDENT = 30        # dies que es continua buscant la imatge de després
VALID_MINIM = 0.85           # part de la zona sense núvols perquè la imatge serveixi
CANDIDATS_MAX = 5            # imatges que es proven per banda
RESOLUCIO = 10               # metres per píxel
DIST_FOCUS_M = 1500          # una taca de crema es compta si toca a menys d'aquesta distància d'un focus

# Llindars de l'índex de crema (dNBR), segons l'escala de l'USGS.
LLINDAR_LLEU = 0.10
LLINDAR_MODERADA = 0.27
LLINDAR_ALTA = 0.66

EVAL_NBR = """//VERSION=3
function setup() {
  return { input: [{ bands: ["B08", "B12", "SCL", "dataMask"] }], output: { bands: 2, sampleType: "FLOAT32" } };
}
function evaluatePixel(s) {
  var nbr = (s.B08 - s.B12) / (s.B08 + s.B12 + 1e-6);
  // Fora: sense dades, saturat, ombra de núvol, núvols, cirrus i neu.
  var dolent = [0, 1, 3, 8, 9, 10, 11].indexOf(s.SCL) !== -1;
  return [nbr, (s.dataMask === 1 && !dolent) ? 1 : 0];
}"""

EVAL_RGB = """//VERSION=3
function setup() { return { input: ["B04", "B03", "B02"], output: { bands: 3 } }; }
function evaluatePixel(s) { return [2.6 * s.B04, 2.6 * s.B03, 2.6 * s.B02]; }"""


# -----------------------------------------------------------------------------
# UTILITATS
# -----------------------------------------------------------------------------
def log(*a):
    print(*a, flush=True)


def data_text(iso):
    d = iso[:10]
    return d[8:10] + "." + d[5:7] + "." + d[0:4]


def tg(metode, dades=None, fitxers=None):
    if not TG_TOKEN:
        log("Telegram: falta el secret TELEGRAM_TOKEN_INFOS.")
        return {}
    try:
        r = requests.post("https://api.telegram.org/bot" + TG_TOKEN + "/" + metode,
                          data=dades or {}, files=fitxers, timeout=120)
        if r.status_code != 200:
            log("Telegram", metode, r.status_code, r.text[:300])
        return r.json()
    except Exception as e:
        log("Telegram", metode, e)
        return {}


_token = {"valor": None, "fins": 0}


def token():
    if _token["valor"] and time.time() < _token["fins"] - 60:
        return _token["valor"]
    if not CLIENT_ID or not CLIENT_SECRET:
        raise RuntimeError("Falten els secrets CDSE_CLIENT_ID i CDSE_CLIENT_SECRET.")
    r = requests.post(URL_TOKEN, data={"grant_type": "client_credentials",
                                       "client_id": CLIENT_ID, "client_secret": CLIENT_SECRET}, timeout=60)
    if r.status_code != 200:
        raise RuntimeError("Copernicus no accepta la clau (HTTP %s): %s" % (r.status_code, r.text[:200]))
    j = r.json()
    _token["valor"] = j["access_token"]
    _token["fins"] = time.time() + int(j.get("expires_in", 600))
    return _token["valor"]


def capcaleres():
    return {"Authorization": "Bearer " + token()}


# -----------------------------------------------------------------------------
# GEOMETRIA (UTM 31N, la de Catalunya)
# -----------------------------------------------------------------------------
def a_utm(lat, lon):
    from pyproj import Transformer
    return Transformer.from_crs("EPSG:4326", "EPSG:32631", always_xy=True).transform(lon, lat)


def zona(lat, lon, radi_km):
    """Quadrat de costat 2*radi al voltant del punt: caixa UTM, caixa geogràfica i mida en píxels."""
    x, y = a_utm(lat, lon)
    m = radi_km * 1000
    caixa_utm = [x - m, y - m, x + m, y + m]
    costat = int(round(2 * m / RESOLUCIO))
    dlat = radi_km / 111.32
    dlon = radi_km / (111.32 * math.cos(math.radians(lat)))
    caixa_geo = [lon - dlon, lat - dlat, lon + dlon, lat + dlat]
    return caixa_utm, caixa_geo, costat


def pixel_de(lat, lon, caixa_utm, costat):
    x, y = a_utm(lat, lon)
    col = (x - caixa_utm[0]) / RESOLUCIO
    fila = (caixa_utm[3] - y) / RESOLUCIO
    return int(fila), int(col)


# -----------------------------------------------------------------------------
# COPERNICUS: CATÀLEG I IMATGES
# -----------------------------------------------------------------------------
def dates_disponibles(caixa_geo, des_de, fins_a):
    """Dates (yyyy-mm-dd) amb imatge Sentinel-2 a la zona, amb núvols de l'escena per sota del 70 %."""
    cos = {"bbox": caixa_geo, "collections": ["sentinel-2-l2a"], "limit": 100,
           "datetime": des_de + "T00:00:00Z/" + fins_a + "T23:59:59Z",
           "filter": "eo:cloud_cover < 70", "filter-lang": "cql2-text",
           "fields": {"include": ["properties.datetime", "properties.eo:cloud_cover"], "exclude": []}}
    r = requests.post(URL_CATALEG, json=cos, headers=capcaleres(), timeout=90)
    if r.status_code != 200:
        raise RuntimeError("Catàleg de Copernicus: HTTP %s %s" % (r.status_code, r.text[:200]))
    dates = {}
    for f in r.json().get("features", []):
        p = f.get("properties", {})
        d = str(p.get("datetime", ""))[:10]
        if d:
            dates[d] = min(dates.get(d, 100), float(p.get("eo:cloud_cover", 100)))
    return dates


def demanar(caixa_utm, costat, dia, evalscript, format_):
    cos = {
        "input": {
            "bounds": {"bbox": caixa_utm, "properties": {"crs": CRS_UTM}},
            "data": [{"type": "sentinel-2-l2a",
                      "dataFilter": {"timeRange": {"from": dia + "T00:00:00Z", "to": dia + "T23:59:59Z"},
                                     "mosaickingOrder": "leastCC"}}]
        },
        "output": {"width": costat, "height": costat,
                   "responses": [{"identifier": "default", "format": {"type": format_}}]},
        "evalscript": evalscript
    }
    r = requests.post(URL_PROCES, json=cos, headers=capcaleres(), timeout=180)
    if r.status_code != 200:
        raise RuntimeError("Copernicus (imatge del %s): HTTP %s %s" % (dia, r.status_code, r.text[:200]))
    return r.content


def nbr_del_dia(caixa_utm, costat, dia):
    import tifffile
    a = tifffile.imread(io.BytesIO(demanar(caixa_utm, costat, dia, EVAL_NBR, "image/tiff")))
    if a.ndim == 3 and a.shape[0] == 2 and a.shape[-1] != 2:
        a = np.moveaxis(a, 0, -1)
    return a[..., 0].astype("float32"), a[..., 1] > 0.5


def millor_imatge(caixa_utm, caixa_geo, costat, des_de, fins_a, ordre):
    """La primera imatge (segons l'ordre) amb prou zona sense núvols. Retorna (dia, nbr, valid, part_valida) o None."""
    dates = dates_disponibles(caixa_geo, des_de, fins_a)
    candidats = sorted(dates, reverse=(ordre == "recent"))[:CANDIDATS_MAX * 2]
    candidats = sorted(candidats, key=lambda d: (dates[d] > 40, candidats.index(d)))[:CANDIDATS_MAX]
    for dia in candidats:
        try:
            nbr, valid = nbr_del_dia(caixa_utm, costat, dia)
        except Exception as e:
            log("  ", dia, e)
            continue
        part = float(valid.mean())
        log("   %s: %d %% de la zona sense núvols" % (dia, round(part * 100)))
        if part >= VALID_MINIM:
            return dia, nbr, valid, part
    return None


# -----------------------------------------------------------------------------
# CÀLCUL DE LA CREMA
# -----------------------------------------------------------------------------
def calcular_crema(nbr_abans, valid_abans, nbr_despres, valid_despres, focus_px):
    from scipy import ndimage
    valid = valid_abans & valid_despres
    dnbr = np.where(valid, nbr_abans - nbr_despres, 0.0)

    # Taques d'afectació (lleu o més) connectades; només compten les que toquen un focus.
    afectat = dnbr >= LLINDAR_LLEU
    etiquetes, n = ndimage.label(afectat, structure=np.ones((3, 3)))
    radi_px = DIST_FOCUS_M / RESOLUCIO
    files, cols = np.indices(dnbr.shape)
    prop_focus = np.zeros(dnbr.shape, bool)
    for (f, c) in focus_px:
        prop_focus |= (files - f) ** 2 + (cols - c) ** 2 <= radi_px ** 2
    bones = set(np.unique(etiquetes[prop_focus & afectat]).tolist()) - {0}
    crema = np.isin(etiquetes, list(bones)) if bones else np.zeros(dnbr.shape, bool)

    ha = RESOLUCIO * RESOLUCIO / 10000.0
    # Perímetre (canvi 05.10.2026): la taca cremada amb les illes de dins
    # omplertes, que és com es calculen les xifres oficials.
    perimetre = ndimage.binary_fill_holes(crema) if crema.any() else crema
    vora = bool(crema[0, :].any() or crema[-1, :].any() or crema[:, 0].any() or crema[:, -1].any())
    lleu = crema & (dnbr < LLINDAR_MODERADA)
    moderada = crema & (dnbr >= LLINDAR_MODERADA) & (dnbr < LLINDAR_ALTA)
    alta = crema & (dnbr >= LLINDAR_ALTA)
    return {
        "dnbr": dnbr, "lleu": lleu, "moderada": moderada, "alta": alta,
        "ha_lleu": float(lleu.sum() * ha), "ha_moderada": float(moderada.sum() * ha),
        "ha_alta": float(alta.sum() * ha), "valid": float(valid.mean()),
        "ha_perimetre": float(perimetre.sum() * ha),
        "ha_sense_dades": float((perimetre & ~valid).sum() * ha),
        "vora": vora
    }


def imatge_rgb(caixa_utm, costat, dia, crema=None, focus_px=None):
    from PIL import Image, ImageDraw
    img = Image.open(io.BytesIO(demanar(caixa_utm, costat, dia, EVAL_RGB, "image/png"))).convert("RGB")
    if crema is not None:
        capa = np.zeros((costat, costat, 4), dtype=np.uint8)
        capa[crema["lleu"]] = (255, 230, 0, 110)
        capa[crema["moderada"]] = (255, 120, 0, 150)
        capa[crema["alta"]] = (220, 0, 0, 180)
        img = Image.alpha_composite(img.convert("RGBA"), Image.fromarray(capa, "RGBA")).convert("RGB")
    if focus_px:
        d = ImageDraw.Draw(img)
        for (f, c) in focus_px:
            d.ellipse([c - 6, f - 6, c + 6, f + 6], outline=(0, 200, 255), width=3)
    if img.width > 1280:
        img = img.resize((1280, 1280))
    sortida = io.BytesIO()
    img.save(sortida, "JPEG", quality=88)
    return sortida.getvalue()


def num(x):
    return ("%.1f" % x).replace(".", ",")


# -----------------------------------------------------------------------------
# INCENDI
# -----------------------------------------------------------------------------
def processar_incendi(inc):
    """Retorna True si s'ha pogut fer (o s'ha de deixar de provar), False si cal tornar-ho a provar."""
    lat, lon, dia = float(inc["lat"]), float(inc["lon"]), inc["data"]
    radi = max(1.0, min(10.0, float(inc.get("radi") or 3)))
    nom = inc.get("nom") or ("Incendi del " + data_text(dia))
    punts = inc.get("punts") or [[lat, lon]]
    log("Incendi:", nom, lat, lon, dia, "radi", radi, "km")

    caixa_utm, caixa_geo, costat = zona(lat, lon, radi)
    focus_px = [pixel_de(p[0], p[1], caixa_utm, costat) for p in punts]
    d_incendi = date.fromisoformat(dia)
    avui = datetime.now(timezone.utc).date()

    log(" Imatge d'abans:")
    abans = millor_imatge(caixa_utm, caixa_geo, costat,
                          (d_incendi - timedelta(days=DIES_ABANS)).isoformat(),
                          (d_incendi - timedelta(days=1)).isoformat(), "recent")
    if not abans:
        tg("sendMessage", {"chat_id": TG_CHAT, "text": "🛰️ " + nom + "\nNo hi ha cap imatge de Sentinel-2 sense núvols "
                           "dels " + str(DIES_ABANS) + " dies anteriors a l'incendi: no es pot calcular la zona cremada."})
        return True

    log(" Imatge de després:")
    if (avui - d_incendi).days < 1:
        return False
    if inc.get("data_fi"):
        # Incendi passat: la primera imatge neta després que s'apagués.
        d_fi = date.fromisoformat(inc["data_fi"])
        despres = millor_imatge(caixa_utm, caixa_geo, costat, (d_fi + timedelta(days=1)).isoformat(),
                                min(avui, d_fi + timedelta(days=60)).isoformat(), "antic")
    else:
        despres = millor_imatge(caixa_utm, caixa_geo, costat, (d_incendi + timedelta(days=1)).isoformat(),
                                avui.isoformat(), "recent")
    if not despres:
        log(" Encara no hi ha cap imatge neta de després.")
        return False

    crema = calcular_crema(abans[1], abans[2], despres[1], despres[2], focus_px)
    total = crema["ha_moderada"] + crema["ha_alta"]
    log(" Crema: lleu %.1f ha · moderada %.1f ha · alta %.1f ha" % (crema["ha_lleu"], crema["ha_moderada"], crema["ha_alta"]))

    img_abans = imatge_rgb(caixa_utm, costat, abans[0], focus_px=focus_px)
    img_despres = imatge_rgb(caixa_utm, costat, despres[0], crema=crema, focus_px=focus_px)

    if total + crema["ha_lleu"] < 1:
        conclusio = ("No s'aprecia cap superfície cremada clara a la zona. Pot ser un foc molt petit, "
                     "sota els arbres, una crema agrícola o una detecció falsa.")
    else:
        illes = max(0.0, crema["ha_perimetre"] - total - crema["ha_lleu"])
        conclusio = ("🔥 Superfície dins del perímetre de l'incendi: " + num(crema["ha_perimetre"]) + " ha\n"
                     "   (comparable amb les xifres oficials, que inclouen les illes sense cremar)\n\n"
                     "Dins del perímetre:\n"
                     "   🔴 gravetat alta: " + num(crema["ha_alta"]) + " ha\n"
                     "   🟠 gravetat moderada: " + num(crema["ha_moderada"]) + " ha\n"
                     "   🟡 afectació lleu: " + num(crema["ha_lleu"]) + " ha\n"
                     "   ⬜ illes sense crema aparent: " + num(illes) + " ha\n"
                     "   ➡️ cremada de manera clara (alta + moderada): " + num(total) + " ha")
        if crema["ha_sense_dades"] >= 1:
            conclusio += "\n☁️ " + num(crema["ha_sense_dades"]) + " ha del perímetre tapades per núvols en alguna de les imatges: la xifra real pot ser més alta."
        if crema["vora"]:
            conclusio += "\n⚠️ La zona cremada arriba a la vora del quadrat analitzat: torna-ho a llançar amb un radi més gran."

    resum = ("🛰️ ZONA CREMADA — " + nom + "\n\n" + conclusio + "\n\n"
             "📷 Imatges de Sentinel-2: abans, " + data_text(abans[0]) + " · després, " + data_text(despres[0]) + "\n"
             "📍 Zona analitzada: " + num(2 * radi) + " × " + num(2 * radi) + " km al voltant de "
             "https://www.google.com/maps?q=%.5f,%.5f\n" % (lat, lon) +
             "ℹ️ Estimació automàtica (índex de crema dNBR): cal contrastar-la. Els camps segats o llaurats entre "
             "les dues dates també poden sortir com a «cremats».\n"
             "🔎 Contrast oficial (incendis de més de 30 ha): https://forest-fire.emergency.copernicus.eu/apps/effis_current_situation/\n"
             "🗺️ Explorar les imatges: https://browser.dataspace.copernicus.eu/?zoom=13&lat=%.5f&lng=%.5f" % (lat, lon))

    media = [{"type": "photo", "media": "attach://abans", "caption": "Abans · " + data_text(abans[0])},
             {"type": "photo", "media": "attach://despres",
              "caption": "Després · " + data_text(despres[0]) + " · zona cremada pintada (vermell: alta, taronja: moderada, groc: lleu)"}]
    tg("sendMediaGroup", {"chat_id": TG_CHAT, "media": json.dumps(media)},
       {"abans": ("abans.jpg", img_abans), "despres": ("despres.jpg", img_despres)})
    tg("sendMessage", {"chat_id": TG_CHAT, "text": resum, "disable_web_page_preview": "true"})
    return True


# -----------------------------------------------------------------------------
# PENDENTS
# -----------------------------------------------------------------------------
def llegir_pendents():
    try:
        with open(PENDENTS, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return []


def desar_pendents(llista):
    with open(PENDENTS, "w", encoding="utf-8") as f:
        json.dump(llista, f, ensure_ascii=False, indent=1)


def revisar_pendents(nous=None):
    llista = llegir_pendents()
    for inc in (nous or []):
        clau = "%.3f|%.3f|%s" % (float(inc["lat"]), float(inc["lon"]), inc["data"])
        if any("%.3f|%.3f|%s" % (float(p["lat"]), float(p["lon"]), p["data"]) == clau for p in llista):
            log("Ja és a la llista de pendents:", clau)
            continue
        inc["afegit"] = datetime.now(timezone.utc).date().isoformat()
        llista.append(inc)

    queden = []
    avui = datetime.now(timezone.utc).date()
    for inc in llista:
        try:
            fet = processar_incendi(inc)
        except Exception as e:
            log("❌", e)
            fet = False
        if fet:
            continue
        if (avui - date.fromisoformat(inc["data"])).days > DIES_MAX_PENDENT:
            tg("sendMessage", {"chat_id": TG_CHAT, "text": "🛰️ " + (inc.get("nom") or "Incendi") + "\nEn " +
                               str(DIES_MAX_PENDENT) + " dies no hi ha hagut cap imatge de Sentinel-2 sense núvols: "
                               "no es pot calcular la zona cremada."})
            continue
        queden.append(inc)
    desar_pendents(queden)
    log("Pendents:", len(queden))


if __name__ == "__main__":
    ordre = sys.argv[1] if len(sys.argv) > 1 else "pendents"
    if ordre == "incendi":
        punts = []
        try:
            punts = json.loads(os.environ.get("PUNTS") or "[]")
        except Exception:
            pass
        inc = {"lat": os.environ["LAT"], "lon": os.environ["LON"],
               "data": (os.environ.get("DATA") or "").strip() or datetime.now(timezone.utc).date().isoformat(),
               "data_fi": (os.environ.get("DATA_FI") or "").strip(),
               "nom": os.environ.get("NOM", "").strip(), "radi": os.environ.get("RADI", "").strip() or "3",
               "punts": punts}
        revisar_pendents([inc])
    else:
        revisar_pendents()
