import pandas as pd
from flask import Flask, request, jsonify
import os
import re

# =========================================================
# CONFIG
# =========================================================

BASE_PATH = os.path.dirname(__file__)
DATA_PATH = os.path.join(BASE_PATH, "data")

PORT = 8000

app = Flask(__name__)

# =========================================================
# OUTILS CIM
# =========================================================

def normalize_cim(code):

    if code is None:
        return None

    code = str(code).strip().upper()

    return code.replace(".", "")


def cim_variants(code):

    code = str(code).strip().upper()

    no_dot = code.replace(".", "")

    if "." in code:
        dot = code
    else:
        if len(no_dot) > 3:
            dot = no_dot[:3] + "." + no_dot[3:]
        else:
            dot = no_dot

    return {no_dot, dot}


def normalize(x):

    if pd.isna(x):
        return None

    return str(x).strip().upper()


def read_csv_safe(path):

    return pd.read_csv(
        path,
        encoding="latin-1",
        sep=None,
        engine="python"
    )

print("Chargement tables PMSI...")

# =========================================================
# CMA
# =========================================================

cma = read_csv_safe(
    os.path.join(DATA_PATH, "tabula-vol1an4.csv")
)

# renommer colonnes tabula

cma = cma.rename(columns={
    "Unnamed: 0": "cim",
    "Unnamed: 1": "niveau",
    "Liste": "liste_diag",
    "Liste.1": "liste_racine",
    "Unnamed: 4": "libelle"
})

# supprimer lignes parasites

cma = cma[cma["niveau"].astype(str).str.contains(r"\d")]

cma["niveau"] = cma["niveau"].astype(int)

# normalisation CIM

cma["cim"] = (
    cma["cim"]
    .astype(str)
    .str.upper()
    .str.replace(".", "", regex=False)
)

# index CMA

cma_dict = {}

for _, row in cma.iterrows():

    cim = row["cim"]

    entry = {
        "niveau": row["niveau"],
        "liste_diag": str(row["liste_diag"]),
        "liste_racine": str(row["liste_racine"])
    }

    if cim not in cma_dict:
        cma_dict[cim] = []

    cma_dict[cim].append(entry)

print("CMA indexées :", len(cma_dict))

# =========================================================
# LISTES EXCLUSION
# =========================================================

exclusion_lists = {}

with open(
    os.path.join(DATA_PATH, "tabula-vol1an5.csv"),
    encoding="latin-1"
) as f:

    for line in f:

        parts = line.strip().split(",")

        if len(parts) < 2:
            continue

        liste = parts[0].strip()

        codes = re.split(r"\s+", " ".join(parts[1:]))

        exclusion_lists[liste] = set(
            [normalize_cim(c) for c in codes if c]
        )

# =========================================================
# CMD
# =========================================================

cmd_df = read_csv_safe(
    os.path.join(DATA_PATH, "tabula-vol1an7.csv")
)

cmd_map = {}

for _, row in cmd_df.iterrows():

    cim = normalize_cim(row.iloc[0])
    cmd = str(row.iloc[1])

    cmd_map[cim] = cmd

# =========================================================
# RACINES
# =========================================================

racines = {}

with open(
    os.path.join(
        DATA_PATH,
        "LISTE DES RACINES DE GHM DE LA CMD n°20.csv"
    ),
    encoding="latin-1"
) as f:

    for line in f:

        parts = line.strip().split(";")

        if len(parts) < 1:
            continue

        racine = parts[0]

        racines[racine] = racine

print("Tables chargées.")

# =========================================================
# FONCTIONS
# =========================================================

def determine_cmd(dp):

    dp = normalize_cim(dp)

    return cmd_map.get(dp)


def determine_racine(dp):

    dp = normalize_cim(dp)

    if dp == "Z502":
        return "20Z04"

    if dp == "Z503":
        return "20Z05"

    return None


def get_cma(code):

    variants = cim_variants(code)

    results = []

    for v in variants:

        key = normalize_cim(v)

        if key in cma_dict:
            results.extend(cma_dict[key])

    return results if results else None


def check_exclusion(dp, cmd, racine, liste_diag, liste_racine):

    excl = False

    dp = normalize_cim(dp)

    if liste_diag in exclusion_lists:

        liste = exclusion_lists[liste_diag]

        if dp in liste or cmd in liste:
            excl = True

    if liste_racine in exclusion_lists:

        liste = exclusion_lists[liste_racine]

        if racine in liste:
            excl = True

    return excl

# =========================================================
# CALCUL CMA
# =========================================================

def compute_cma(dp, cmd, racine, codes):

    niveaux = []
    details = []

    for code in codes:

        entries = get_cma(code)

        if not entries:

            details.append({
                "code": code,
                "status": "NON_CMA"
            })

            continue

        active = False
        niveau_max = 0

        for entry in entries:

            excl = check_exclusion(
                dp,
                cmd,
                racine,
                entry["liste_diag"],
                entry["liste_racine"]
            )

            if not excl:

                active = True
                niveau_max = max(
                    niveau_max,
                    entry["niveau"]
                )

        if active:

            niveaux.append(niveau_max)

            details.append({
                "code": code,
                "status": "CMA_ACTIVE",
                "niveau": niveau_max
            })

        else:

            details.append({
                "code": code,
                "status": "CMA_NEUTRALISEE"
            })

    niveau_final = max(niveaux) if niveaux else 1

    return niveau_final, details

# =========================================================
# CMA FREQUEMMENT ASSOCIEES (LOGIQUE CLINIQUE)
# =========================================================

CLINICAL_CMA_PATTERNS = {

    # cirrhose
    "K703": ["K766","R18","D696","K729"],

    # varices oesophagiennes
    "I850": ["K766"],

    # alcool
    "F106": ["G312","G621"],

}

# =========================================================
# SUGGESTION CMA POTENTIELLES
# =========================================================

def suggest_potential_cma(codes):

    suggestions = set()

    normalized = [normalize_cim(c) for c in codes]

    for code in normalized:

        if code in CLINICAL_CMA_PATTERNS:

            for cma in CLINICAL_CMA_PATTERNS[code]:

                if cma not in normalized:

                    suggestions.add(cma)

    return list(suggestions)

# =========================================================
# EXPLORATION DP ALTERNATIFS
# =========================================================

DP_ALTERNATIVES = {

    "Z502": ["F103","K703"],
    "Z503": ["F113","F192"]

}


def explore_dp(dp, codes):

    results = []

    normalized_dp = normalize_cim(dp)

    if normalized_dp not in DP_ALTERNATIVES:
        return results

    for alt_dp in DP_ALTERNATIVES[normalized_dp]:

        cmd = determine_cmd(alt_dp)

        racine = determine_racine(alt_dp)

        if racine is None:
            continue

        niveau, _ = compute_cma(
            alt_dp,
            cmd,
            racine,
            codes
        )

        ghm = racine + str(niveau)

        results.append({

            "dp": alt_dp,
            "cmd": cmd,
            "racine": racine,
            "niveau": niveau,
            "ghm": ghm

        })

    return results

# =========================================================
# GROUPEUR AVANCE
# =========================================================

def grouper(dp, codes):

    cmd = determine_cmd(dp)

    racine = determine_racine(dp)

    niveau, details = compute_cma(
        dp,
        cmd,
        racine,
        codes
    )

    ghm = racine + str(niveau)

    # suggestions cliniques
    suggestions = explore_clinical_codes(codes)

    # CMA potentielles
    potential_cma = suggest_potential_cma(codes)

    # exploration DP
    dp_exploration = explore_dp(dp, codes)

    return {

        "dp": dp,
        "cmd": cmd,
        "racine": racine,
        "niveau": niveau,
        "ghm": ghm,

        "details": details,

        "clinical_suggestions": suggestions,

        "potential_cma": potential_cma,

        "dp_exploration": dp_exploration

    }

# =========================================================
# EXPLORATION CIM CLINIQUEMENT PROCHES
# =========================================================

CLINICAL_EQUIVALENTS = {

    # varices oesophagiennes
    "I859": ["I850"],

    # cirrhose
    "K703": ["K746", "K704"],

    # trouble cognitif alcool
    "F067": ["F106"],

    # dépression sévère
    "F322": ["F332"]

}


def explore_clinical_codes(codes):

    suggestions = {}

    for code in codes:

        key = normalize_cim(code)

        if key in CLINICAL_EQUIVALENTS:

            suggestions[code] = CLINICAL_EQUIVALENTS[key]

    return suggestions


# =========================================================
# API
# =========================================================

@app.route("/groupage", methods=["POST"])
def api_groupage():

    data = request.json

    dp = data["dp"]
    codes = data["codes"]

    result = grouper(dp, codes)

    return jsonify(result)


# =========================================================
# TEST
# =========================================================

@app.route("/test")
def test():

    test_case = {

        "dp": "Z50.2",
        "codes": [
            "K70.3",
            "I85.9",
            "F32.2",
            "F06.7"
        ]

    }

    result = grouper(
        test_case["dp"],
        test_case["codes"]
    )

    return jsonify(result)


# =========================================================
# START SERVEUR
# =========================================================

if __name__ == "__main__":

    print("Moteur PMSI prêt")

    app.run(
        host="0.0.0.0",
        port=PORT,
        debug=True
    )