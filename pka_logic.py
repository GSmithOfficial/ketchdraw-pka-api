"""
pKa calculation logic using UnipKa.
Apache-2.0 License - uses finlayiainmaclean/unipka wrapper.
"""
from unipka import UnipKa
from rdkit import Chem
import math

# Initialize model once (loaded when module imports)
print("🔄 Loading UnipKa model...")
calc = UnipKa()
print("✅ UnipKa model loaded")


def calculate_pka_properties(smiles: str, pH: float = 7.4) -> dict:
    """
    Calculate pKa and related properties for a molecule.
    Returns a JSON-serializable dictionary.
    """
    result = {
        "input_smiles": smiles,
        "query_pH": pH,
        "success": False,
        "error": None,
        "pka": {},
        "logd": None,
        "state_penalty": None,
        "dominant_microstate": None,
        "microstates": []
    }
    
    try:
        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["error"] = "Invalid SMILES string"
            return result
        
        # Canonical SMILES
        canonical_smiles = Chem.MolToSmiles(mol)
        result["canonical_smiles"] = canonical_smiles
        
        # pKa values (handle nan)
        acidic = calc.get_acidic_macro_pka(smiles)
        basic = calc.get_basic_macro_pka(smiles)
        result["pka"]["acidic"] = None if math.isnan(acidic) else round(acidic, 2)
        result["pka"]["basic"] = None if math.isnan(basic) else round(basic, 2)
        
        # logD at query pH
        logd = calc.get_logd(smiles, pH=pH)
        result["logd"] = round(logd, 3)
        
        # State penalty
        penalty_tuple = calc.get_state_penalty(smiles, pH=pH)
        result["state_penalty"] = round(penalty_tuple[0], 6)
        
        # Dominant microstate at query pH
        dominant = calc.get_dominant_microstate(smiles, pH=pH)
        result["dominant_microstate"] = Chem.MolToSmiles(dominant)
        
        # Full distribution
        dist_df = calc.get_distribution(smiles, pH=pH)
        for _, row in dist_df.iterrows():
            result["microstates"].append({
                "charge": int(row["charge"]),
                "smiles": row["smiles"],
                "population": round(float(row["population"]), 8),
                "free_energy": round(float(row["free_energy"]), 3),
                "is_query_mol": bool(row["is_query_mol"])
            })
        
        result["success"] = True
        
    except Exception as e:
        result["error"] = str(e)
    
    return result
