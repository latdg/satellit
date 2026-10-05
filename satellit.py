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
# MÒDUL 2: ZONES INUNDADES AMB SENTINEL-1 (radar) — afegit 05.10.2026
#   Per a un episodi de pluja (data), agafa la primera passada del radar el mateix
#   dia o els dies següents i la compara amb una passada anterior en sec (de la
#   mateixa òrbita). L'aigua estesa torna molt poc senyal: els píxels que eren
#   secs i ara retornen molt poc senyal es compten com a inundats. El radar veu
#   igual amb núvols i de nit. Si en INUNDACIO_DIES_MAX dies no hi ha cap passada,
#   s'avisa i es deixa córrer (l'aigua ja hauria baixat).
#
# MÒDULS 3 I 4: INFORMES MENSUALS AMB SENTINEL-2 — afegits 05.10.2026
#   Cada mes (dia 5) es fa una imatge "neta" del mes anterior: per a cada píxel,
#   la mitjana (mediana) de totes les passades sense núvols del mes.
#   3. VEGETACIÓ: l'índex de verdor (NDVI) del mes comparat amb el mateix mes dels
#      anys anteriors (des del 2019): on la vegetació està més seca o més verda
#      que de costum (sequera, risc d'incendi, recuperació després de pluges).
#   4. CANVIS AL TERRITORI: zones que fa un any tenien vegetació i ara no (tales,
#      obres, pedreres, urbanitzacions, cremes, també collites), amb les imatges
#      d'abans i d'ara de les més grans, per contrastar-les amb llicències, tauler
#      d'anuncis i contractació.
#
# Ordres:
#   python3 satellit.py mensual    (vegetació i canvis del mes anterior; DATA=aaaa-mm per a un altre mes)
#   python3 satellit.py vegetacio  /  python3 satellit.py canvis
#   python3 satellit.py inundacio  (dades a les variables DATA, NOM; LAT, LON i RADI opcionals)
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
# INUNDACIONS (SENTINEL-1)
# -----------------------------------------------------------------------------
INUNDACIO_CENTRE = (41.9597, 3.0386)      # la Bisbal d'Empordà
INUNDACIO_RADI_KM = 12                    # quadrat de 24 x 24 km (els quatre municipis i la plana del Daró)
INUNDACIO_RESOLUCIO = 20                  # metres per píxel (prou per a camps negats)
INUNDACIO_DIES_MAX = 4                    # dies després de la pluja que es busca una passada
AIGUA_DB = -18.0                          # per sota d'aquest senyal (dB), aigua
CAIGUDA_DB = 3.0                          # el senyal ha de baixar com a mínim això respecte de la referència
MIN_PIXELS_TACA = 5                       # taques més petites, soroll

EVAL_VV = """//VERSION=3
function setup() { return { input: [{ bands: ["VV", "dataMask"] }], output: { bands: 2, sampleType: "FLOAT32" } }; }
function evaluatePixel(s) { return [s.VV, s.dataMask]; }"""


def passades_s1(caixa_geo, des_de, fins_a):
    """[(dia, orbita)] de les passades de Sentinel-1 a la zona."""
    cos = {"bbox": caixa_geo, "collections": ["sentinel-1-grd"], "limit": 100,
           "datetime": des_de + "T00:00:00Z/" + fins_a + "T23:59:59Z",
           "fields": {"include": ["properties.datetime", "properties.sat:orbit_state"], "exclude": []}}
    r = requests.post(URL_CATALEG, json=cos, headers=capcaleres(), timeout=90)
    if r.status_code != 200:
        raise RuntimeError("Catàleg de Copernicus (radar): HTTP %s %s" % (r.status_code, r.text[:200]))
    vistes = {}
    for f in r.json().get("features", []):
        p = f.get("properties", {})
        dia = str(p.get("datetime", ""))[:10]
        orbita = str(p.get("sat:orbit_state", "")).upper()
        if dia and orbita in ("ASCENDING", "DESCENDING"):
            vistes[(dia, orbita)] = True
    return sorted(vistes)


def vv_db(caixa_utm, costat, dia, orbita):
    import tifffile
    from scipy import ndimage
    cos = {
        "input": {
            "bounds": {"bbox": caixa_utm, "properties": {"crs": CRS_UTM}},
            "data": [{"type": "sentinel-1-grd",
                      "dataFilter": {"timeRange": {"from": dia + "T00:00:00Z", "to": dia + "T23:59:59Z"},
                                     "acquisitionMode": "IW", "polarization": "DV", "orbitDirection": orbita},
                      "processing": {"backCoeff": "GAMMA0_TERRAIN", "orthorectify": True, "demInstance": "COPERNICUS"}}]
        },
        "output": {"width": costat, "height": costat,
                   "responses": [{"identifier": "default", "format": {"type": "image/tiff"}}]},
        "evalscript": EVAL_VV
    }
    r = requests.post(URL_PROCES, json=cos, headers=capcaleres(), timeout=180)
    if r.status_code != 200:
        raise RuntimeError("Copernicus (radar del %s): HTTP %s %s" % (dia, r.status_code, r.text[:200]))
    a = tifffile.imread(io.BytesIO(r.content))
    if a.ndim == 3 and a.shape[0] == 2 and a.shape[-1] != 2:
        a = np.moveaxis(a, 0, -1)
    vv, valid = a[..., 0], a[..., 1] > 0.5
    db = 10 * np.log10(np.clip(vv, 1e-5, None))
    db = ndimage.median_filter(db, size=3)        # treu el "soroll" granulat del radar
    return db.astype("float32"), valid


def calcular_inundacio(ref_db, ref_valid, post_db, post_valid):
    from scipy import ndimage
    valid = ref_valid & post_valid
    aigua_abans = ref_db < AIGUA_DB                   # rius, basses i mar: ja hi eren
    inundat = valid & (post_db < AIGUA_DB) & ((ref_db - post_db) >= CAIGUDA_DB) & ~aigua_abans
    etiquetes, n = ndimage.label(inundat, structure=np.ones((3, 3)))
    if n:
        mides = ndimage.sum(inundat, etiquetes, index=np.arange(1, n + 1))
        grans = np.arange(1, n + 1)[mides >= MIN_PIXELS_TACA]
        inundat = np.isin(etiquetes, grans)
        etiquetes, n = ndimage.label(inundat, structure=np.ones((3, 3)))
    ha_pixel = INUNDACIO_RESOLUCIO * INUNDACIO_RESOLUCIO / 10000.0
    taques = []
    if n:
        mides = ndimage.sum(inundat, etiquetes, index=np.arange(1, n + 1))
        centres = ndimage.center_of_mass(inundat, etiquetes, index=np.arange(1, n + 1))
        for i in np.argsort(mides)[::-1][:5]:
            taques.append({"ha": float(mides[i] * ha_pixel), "fila": centres[i][0], "col": centres[i][1]})
    return {"mascara": inundat, "ha": float(inundat.sum() * ha_pixel), "taques": taques,
            "valid": float(valid.mean())}


def imatge_radar(db, mascara=None):
    from PIL import Image
    gris = np.clip((db + 25) / 25 * 255, 0, 255).astype(np.uint8)
    img = Image.fromarray(gris, "L").convert("RGB")
    if mascara is not None:
        capa = np.zeros(db.shape + (4,), dtype=np.uint8)
        capa[mascara] = (0, 120, 255, 200)
        img = Image.alpha_composite(img.convert("RGBA"), Image.fromarray(capa, "RGBA")).convert("RGB")
    if img.width > 1280:
        img = img.resize((1280, 1280))
    sortida = io.BytesIO()
    img.save(sortida, "JPEG", quality=88)
    return sortida.getvalue()


def lat_lon_de(fila, col, caixa_utm, resolucio):
    from pyproj import Transformer
    x = caixa_utm[0] + (col + 0.5) * resolucio
    y = caixa_utm[3] - (fila + 0.5) * resolucio
    lon, lat = Transformer.from_crs("EPSG:32631", "EPSG:4326", always_xy=True).transform(x, y)
    return lat, lon


def processar_inundacio(inc):
    lat = float(inc.get("lat") or INUNDACIO_CENTRE[0])
    lon = float(inc.get("lon") or INUNDACIO_CENTRE[1])
    radi = max(3.0, min(12.0, float(inc.get("radi") or INUNDACIO_RADI_KM)))
    dia = inc["data"]
    nom = inc.get("nom") or ("Pluges del " + data_text(dia))
    log("Inundació:", nom, dia)

    x, y = a_utm(lat, lon)
    m = radi * 1000
    caixa_utm = [x - m, y - m, x + m, y + m]
    costat = int(round(2 * m / INUNDACIO_RESOLUCIO))
    dlat, dlon = radi / 111.32, radi / (111.32 * math.cos(math.radians(lat)))
    caixa_geo = [lon - dlon, lat - dlat, lon + dlon, lat + dlat]

    d_pluja = date.fromisoformat(dia)
    avui = datetime.now(timezone.utc).date()
    limit = min(avui, d_pluja + timedelta(days=INUNDACIO_DIES_MAX))

    despres = passades_s1(caixa_geo, dia, limit.isoformat())
    if not despres:
        if (avui - d_pluja).days > INUNDACIO_DIES_MAX:
            tg("sendMessage", {"chat_id": TG_CHAT, "text": "🛰️ " + nom + "\nEl radar Sentinel-1 no ha passat per la zona "
                               "en els " + str(INUNDACIO_DIES_MAX) + " dies posteriors a la pluja: no hi ha mapa d'inundació."})
            return True
        log(" Encara no hi ha cap passada del radar després de la pluja.")
        return False
    dia_post, orbita = despres[0]

    # Referència en sec: la passada més recent de la mateixa òrbita, entre 5 i 40 dies abans.
    abans = [p for p in passades_s1(caixa_geo, (d_pluja - timedelta(days=40)).isoformat(),
                                    (d_pluja - timedelta(days=5)).isoformat()) if p[1] == orbita]
    if not abans:
        tg("sendMessage", {"chat_id": TG_CHAT, "text": "🛰️ " + nom + "\nNo hi ha cap passada anterior del radar de la mateixa "
                           "òrbita per comparar: no es pot fer el mapa d'inundació."})
        return True
    dia_ref = abans[-1][0]
    log(" Radar: referència", dia_ref, "· després", dia_post, "·", orbita)

    ref_db, ref_valid = vv_db(caixa_utm, costat, dia_ref, orbita)
    post_db, post_valid = vv_db(caixa_utm, costat, dia_post, orbita)
    res = calcular_inundacio(ref_db, ref_valid, post_db, post_valid)
    log(" Inundat: %.1f ha" % res["ha"])

    linies = ["🛰️ ZONES INUNDADES (radar Sentinel-1) — " + nom, ""]
    if res["ha"] < 1:
        linies.append("💧 El radar no hi veu cap zona negada clara a l'hora de la passada (" + data_text(dia_post) + ").")
    else:
        linies.append("💧 Superfície amb aigua estesa: " + num(res["ha"]) + " ha")
        linies.append("Zones més grans:")
        for t in res["taques"]:
            la, lo = lat_lon_de(t["fila"], t["col"], caixa_utm, INUNDACIO_RESOLUCIO)
            linies.append("   • " + num(t["ha"]) + " ha · https://www.google.com/maps?q=%.5f,%.5f" % (la, lo))
    linies += ["",
               "📷 Radar del " + data_text(dia_post) + " comparat amb el del " + data_text(dia_ref) + " (en sec, mateixa òrbita)",
               "📍 Zona analitzada: " + num(2 * radi) + " × " + num(2 * radi) + " km al voltant de la Bisbal",
               "ℹ️ És la situació a l'hora de la passada: si l'aigua ja havia baixat, no surt. El radar no veu bé "
               "l'aigua als carrers ni sota els arbres. Estimació automàtica: cal contrastar-la.",
               "🗺️ Explorar: https://browser.dataspace.copernicus.eu/?zoom=12&lat=%.5f&lng=%.5f" % (lat, lon)]

    media = [{"type": "photo", "media": "attach://ref", "caption": "Radar en sec · " + data_text(dia_ref)},
             {"type": "photo", "media": "attach://post",
              "caption": "Radar del " + data_text(dia_post) + " · en blau, l'aigua nova"}]
    tg("sendMediaGroup", {"chat_id": TG_CHAT, "media": json.dumps(media)},
       {"ref": ("ref.jpg", imatge_radar(ref_db)), "post": ("post.jpg", imatge_radar(post_db, res["mascara"]))})
    tg("sendMessage", {"chat_id": TG_CHAT, "text": "\n".join(linies), "disable_web_page_preview": "true"})
    return True


# -----------------------------------------------------------------------------
# INFORMES MENSUALS (SENTINEL-2): VEGETACIÓ I CANVIS AL TERRITORI
# -----------------------------------------------------------------------------
MENSUAL_RADI_KM = 12            # 24 x 24 km: els quatre municipis i bona part de les Gavarres
MENSUAL_RESOLUCIO = 20          # metres per píxel
ANY_INICI_HISTORIC = 2019       # primer any de referència per a la vegetació
CANVI_NDVI = 0.30               # caiguda de verdor que es compta com a canvi
CANVI_NDVI_ABANS = 0.40         # només zones que abans tenien vegetació clara
CANVI_MIN_HA = 0.5              # taques més petites no s'avisen
CANVIS_MAX = 8                  # zones que es llisten
MESOS_CAT = ["", "gener", "febrer", "març", "abril", "maig", "juny", "juliol", "agost",
             "setembre", "octubre", "novembre", "desembre"]

EVAL_COMPOSICIO = """//VERSION=3
function setup() {
  return { input: [{ bands: ["B04", "B08", "B11", "SCL", "dataMask"] }],
           output: { bands: 3, sampleType: "FLOAT32" }, mosaicking: "ORBIT" };
}
function mitjana(a) {
  if (!a.length) return NaN;
  a.sort(function (x, y) { return x - y; });
  var m = Math.floor(a.length / 2);
  return a.length % 2 ? a[m] : (a[m - 1] + a[m]) / 2;
}
function evaluatePixel(mostres) {
  var n = [], w = [];
  for (var i = 0; i < mostres.length; i++) {
    var x = mostres[i];
    if (x.dataMask !== 1 || [0, 1, 3, 8, 9, 10, 11].indexOf(x.SCL) !== -1) continue;
    n.push((x.B08 - x.B04) / (x.B08 + x.B04 + 1e-6));
    w.push(x.B11);
  }
  return [mitjana(n), mitjana(w), n.length];
}"""


def mes_anterior():
    avui = datetime.now(timezone.utc).date().replace(day=1)
    d = avui - timedelta(days=1)
    return d.year, d.month


def limits_mes(any_, mes):
    ini = date(any_, mes, 1)
    fi = (date(any_ + (mes == 12), mes % 12 + 1, 1) - timedelta(days=1))
    return ini.isoformat(), fi.isoformat()


def zona_mensual():
    x, y = a_utm(INUNDACIO_CENTRE[0], INUNDACIO_CENTRE[1])
    m = MENSUAL_RADI_KM * 1000
    return [x - m, y - m, x + m, y + m], int(round(2 * m / MENSUAL_RESOLUCIO))


def demanar_rang(caixa_utm, costat, des, fins, evalscript, format_, ordre="leastCC"):
    cos = {
        "input": {
            "bounds": {"bbox": caixa_utm, "properties": {"crs": CRS_UTM}},
            "data": [{"type": "sentinel-2-l2a",
                      "dataFilter": {"timeRange": {"from": des + "T00:00:00Z", "to": fins + "T23:59:59Z"},
                                     "maxCloudCoverage": 80, "mosaickingOrder": ordre}}]
        },
        "output": {"width": costat, "height": costat,
                   "responses": [{"identifier": "default", "format": {"type": format_}}]},
        "evalscript": evalscript
    }
    r = requests.post(URL_PROCES, json=cos, headers=capcaleres(), timeout=300)
    if r.status_code != 200:
        raise RuntimeError("Copernicus (" + des + " - " + fins + "): HTTP %s %s" % (r.status_code, r.text[:200]))
    return r.content


def composicio(caixa_utm, costat, any_, mes):
    """Imatge neta del mes: [verdor, SWIR, passades vàlides] per píxel (NaN si no n'hi ha cap)."""
    import tifffile
    des, fins = limits_mes(any_, mes)
    a = tifffile.imread(io.BytesIO(demanar_rang(caixa_utm, costat, des, fins, EVAL_COMPOSICIO, "image/tiff")))
    if a.ndim == 3 and a.shape[0] == 3 and a.shape[-1] != 3:
        a = np.moveaxis(a, 0, -1)
    return a[..., 0].astype("float32"), a[..., 1].astype("float32")


def num_pct(x):
    return ("%+.0f" % x).replace("-", "−") + " %"


def mapa_anomalia(anom):
    from PIL import Image
    img = np.full(anom.shape + (3,), 235, dtype=np.uint8)
    valid = ~np.isnan(anom)
    a = np.clip(np.nan_to_num(anom) / 0.25, -1, 1)
    sec = valid & (a < 0)
    verd = valid & (a >= 0)
    img[sec, 0] = 255
    img[sec, 1] = (255 * (1 + a[sec])).astype(np.uint8)
    img[sec, 2] = (255 * (1 + a[sec])).astype(np.uint8)
    img[verd, 0] = (255 * (1 - a[verd])).astype(np.uint8)
    img[verd, 1] = 255 - (80 * a[verd]).astype(np.uint8)
    img[verd, 2] = (255 * (1 - a[verd])).astype(np.uint8)
    pil = Image.fromarray(img, "RGB")
    if pil.width > 1200:
        pil = pil.resize((1200, 1200))
    sortida = io.BytesIO()
    pil.save(sortida, "JPEG", quality=88)
    return sortida.getvalue()


def informe_vegetacio(any_, mes):
    caixa, costat = zona_mensual()
    nom_mes = MESOS_CAT[mes] + " de " + str(any_)
    log("Vegetació:", nom_mes)
    ndvi, _ = composicio(caixa, costat, any_, mes)

    historic = []
    for a in range(ANY_INICI_HISTORIC, any_):
        try:
            historic.append(composicio(caixa, costat, a, mes)[0])
        except Exception as e:
            log("  ", a, e)
    if not historic or np.all(np.isnan(ndvi)):
        tg("sendMessage", {"chat_id": TG_CHAT, "text": "🌿 Estat de la vegetació — " + nom_mes +
                           "\nNo hi ha prou imatges sense núvols per fer l'informe d'aquest mes."})
        return
    ref = np.nanmedian(np.stack(historic), axis=0)
    anom = ndvi - ref
    valid = ~np.isnan(anom)
    bosc = valid & (ref >= 0.6)

    def mitjana(m):
        return float(np.nanmean(ndvi[m])), float(np.nanmean(ref[m]))

    tot_ara, tot_ref = mitjana(valid)
    linies = ["🌿 ESTAT DE LA VEGETACIÓ — " + nom_mes,
              "Zona: 24 × 24 km al voltant de la Bisbal (els quatre municipis i bona part de les Gavarres)", "",
              "Índex de verdor (NDVI) mitjà: " + ("%.2f" % tot_ara).replace(".", ",") +
              " (" + num_pct((tot_ara - tot_ref) / tot_ref * 100) + " respecte a la mitjana de " + MESOS_CAT[mes] +
              " de " + str(ANY_INICI_HISTORIC) + "-" + str(any_ - 1) + ")"]
    if bosc.sum() > 100:
        b_ara, b_ref = mitjana(bosc)
        linies.append("🌲 Boscos (zones amb vegetació densa): " + num_pct((b_ara - b_ref) / b_ref * 100))
    sec = float((valid & (anom <= -0.1)).sum() / max(1, valid.sum()) * 100)
    verd = float((valid & (anom >= 0.1)).sum() / max(1, valid.sum()) * 100)
    linies += ["🟥 Superfície clarament més seca que de costum: " + ("%.0f" % sec) + " %",
               "🟩 Superfície clarament més verda que de costum: " + ("%.0f" % verd) + " %"]
    if np.isnan(ndvi).mean() > 0.2:
        linies.append("☁️ Un " + ("%.0f" % (np.isnan(ndvi).mean() * 100)) + " % de la zona no té cap imatge sense núvols aquest mes.")
    linies += ["",
               "ℹ️ Un índex més baix pot indicar sequera, però també incendis, tales o collites. "
               "Mapa: vermell = més sec que de costum; verd = més verd. Estimació automàtica: cal contrastar-la.",
               "Font: Sentinel-2 (Copernicus)"]
    tg("sendPhoto", {"chat_id": TG_CHAT, "caption": "🌿 Vegetació — " + nom_mes + " (vermell: més sec · verd: més verd)"},
       {"photo": ("vegetacio.jpg", mapa_anomalia(anom))})
    tg("sendMessage", {"chat_id": TG_CHAT, "text": "\n".join(linies), "disable_web_page_preview": "true"})


def informe_canvis(any_, mes):
    from scipy import ndimage
    caixa, costat = zona_mensual()
    nom_mes = MESOS_CAT[mes] + " de " + str(any_)
    nom_abans = MESOS_CAT[mes] + " de " + str(any_ - 1)
    log("Canvis al territori:", nom_mes, "respecte a", nom_abans)
    ndvi, swir = composicio(caixa, costat, any_, mes)
    ndvi0, swir0 = composicio(caixa, costat, any_ - 1, mes)

    valid = ~np.isnan(ndvi) & ~np.isnan(ndvi0)
    perdua = valid & (ndvi0 >= CANVI_NDVI_ABANS) & ((ndvi0 - ndvi) >= CANVI_NDVI)
    etiquetes, n = ndimage.label(perdua, structure=np.ones((3, 3)))
    ha_px = MENSUAL_RESOLUCIO * MENSUAL_RESOLUCIO / 10000.0
    taques = []
    if n:
        idx = np.arange(1, n + 1)
        mides = ndimage.sum(perdua, etiquetes, idx)
        centres = ndimage.center_of_mass(perdua, etiquetes, idx)
        dswir = ndimage.mean(np.nan_to_num(swir - swir0), etiquetes, idx)
        for i in np.argsort(mides)[::-1]:
            ha = float(mides[i] * ha_px)
            if ha < CANVI_MIN_HA or len(taques) >= CANVIS_MAX:
                break
            la, lo = lat_lon_de(centres[i][0], centres[i][1], caixa, MENSUAL_RESOLUCIO)
            taques.append({"ha": ha, "lat": la, "lon": lo,
                           "tipus": "sòl nu o construcció" if dswir[i] >= 0.05 else "pèrdua de vegetació"})

    linies = ["🏗️ CANVIS AL TERRITORI — " + nom_mes + " respecte a " + nom_abans,
              "Zona: 24 × 24 km al voltant de la Bisbal", ""]
    if not taques:
        linies.append("No s'hi aprecia cap canvi gran (de " + num(CANVI_MIN_HA) + " ha o més).")
    else:
        linies.append("Zones que fa un any tenien vegetació i ara no (de més gran a més petita):")
        for k, t in enumerate(taques, 1):
            linies.append("%d. %s ha · %s · https://www.google.com/maps?q=%.5f,%.5f" % (k, num(t["ha"]), t["tipus"], t["lat"], t["lon"]))
        linies += ["",
                   "ℹ️ Poden ser obres, tales, pedreres o urbanitzacions, però també collites, llaurades o cremes. "
                   "Val la pena contrastar-ho amb llicències, el tauler d'anuncis i la contractació.",
                   "Font: Sentinel-2 (Copernicus) · estimació automàtica"]

        # Imatges d'abans i d'ara de les tres zones més grans (1 x 1 km).
        media, fitxers = [], {}
        da, fa = limits_mes(any_ - 1, mes)
        db, fb = limits_mes(any_, mes)
        for k, t in enumerate(taques[:3], 1):
            x, y = a_utm(t["lat"], t["lon"])
            c = [x - 500, y - 500, x + 500, y + 500]
            try:
                fitxers["a%d" % k] = ("a%d.png" % k, demanar_rang(c, 100, da, fa, EVAL_RGB, "image/png"))
                fitxers["b%d" % k] = ("b%d.png" % k, demanar_rang(c, 100, db, fb, EVAL_RGB, "image/png"))
                media.append({"type": "photo", "media": "attach://a%d" % k, "caption": "Zona %d · %s" % (k, nom_abans)})
                media.append({"type": "photo", "media": "attach://b%d" % k, "caption": "Zona %d · %s" % (k, nom_mes)})
            except Exception as e:
                log("  imatge zona", k, e)
        if media:
            tg("sendMediaGroup", {"chat_id": TG_CHAT, "media": json.dumps(media)}, fitxers)
    tg("sendMessage", {"chat_id": TG_CHAT, "text": "\n".join(linies), "disable_web_page_preview": "true"})


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
        clau_de = lambda p: (p.get("tipus") or "incendi") + "|" + ("%.3f|%.3f" % (float(p.get("lat") or 0), float(p.get("lon") or 0))) + "|" + p["data"]
        clau = clau_de(inc)
        if any(clau_de(p) == clau for p in llista):
            log("Ja és a la llista de pendents:", clau)
            continue
        inc["afegit"] = datetime.now(timezone.utc).date().isoformat()
        llista.append(inc)

    queden = []
    avui = datetime.now(timezone.utc).date()
    for inc in llista:
        try:
            fet = processar_inundacio(inc) if inc.get("tipus") == "inundacio" else processar_incendi(inc)
        except Exception as e:
            log("❌", e)
            fet = False
        if fet:
            continue
        if inc.get("tipus") != "inundacio" and (avui - date.fromisoformat(inc["data"])).days > DIES_MAX_PENDENT:
            tg("sendMessage", {"chat_id": TG_CHAT, "text": "🛰️ " + (inc.get("nom") or "Incendi") + "\nEn " +
                               str(DIES_MAX_PENDENT) + " dies no hi ha hagut cap imatge de Sentinel-2 sense núvols: "
                               "no es pot calcular la zona cremada."})
            continue
        queden.append(inc)
    desar_pendents(queden)
    log("Pendents:", len(queden))


if __name__ == "__main__":
    ordre = sys.argv[1] if len(sys.argv) > 1 else "pendents"
    if ordre in ("mensual", "vegetacio", "canvis"):
        d = (os.environ.get("DATA") or "").strip()
        any_, mes = (int(d[:4]), int(d[5:7])) if len(d) >= 7 else mes_anterior()
        for nom_, funcio in (("vegetacio", informe_vegetacio), ("canvis", informe_canvis)):
            if ordre in ("mensual", nom_):
                try:
                    funcio(any_, mes)
                except Exception as e:
                    log("❌", nom_, e)
                    tg("sendMessage", {"chat_id": TG_CHAT, "text": "🛰️ No s'ha pogut fer l'informe de " +
                                       ("vegetació" if nom_ == "vegetacio" else "canvis al territori") + ": " + str(e)[:300]})
    elif ordre == "inundacio":
        revisar_pendents([{"tipus": "inundacio",
                           "data": (os.environ.get("DATA") or "").strip() or datetime.now(timezone.utc).date().isoformat(),
                           "nom": os.environ.get("NOM", "").strip(), "lat": os.environ.get("LAT", "").strip(),
                           "lon": os.environ.get("LON", "").strip(), "radi": os.environ.get("RADI", "").strip()}])
    elif ordre == "incendi":
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
