import pandas as pd
import swifter
from typing import List, Dict, Set


def create_trait_modules(trait_df: pd.DataFrame, 
                        individual_df: pd.DataFrame, 
                        icd10_column: str = 'icd10_codes',
                        trait_index_col: str = 'idx',
                        trait_codes_col: str = 'icd10') -> pd.DataFrame:
    """
    Create trait modules for individuals based on their ICD10 codes.
    
    Parameters:
    -----------
    trait_df : pd.DataFrame
        DataFrame with trait definitions. Should have columns:
        - trait_index_col: trait identifier (e.g., 'idx')
        - trait_codes_col: list of ICD10 codes for each trait
    individual_df : pd.DataFrame
        DataFrame with individuals. Should have:
        - icd10_column: column containing ICD10 codes for each individual
    icd10_column : str
        Name of column in individual_df containing ICD10 codes
    trait_index_col : str
        Name of column in trait_df containing trait identifiers
    trait_codes_col : str
        Name of column in trait_df containing ICD10 code lists
        
    Returns:
    --------
    pd.DataFrame
        Original individual_df with added 'trait_modules' column containing
        list of trait indices for each individual
    """
    
    # Create a copy to avoid modifying original
    result_df = individual_df.copy()
    
    # Create mapping from ICD10 code to trait indices
    icd10_to_traits: Dict[str, List[int]] = {}
    
    for _, row in trait_df.iterrows():
        trait_idx = row[trait_index_col]
        icd10_codes = row[trait_codes_col]
        
        # Handle case where icd10_codes might be a string representation of list
        if isinstance(icd10_codes, str):
            # Parse string like "[I48, I489, I482]" to list
            icd10_codes = icd10_codes.strip('[]').replace(' ', '').split(',')
        
        # Add this trait to each ICD10 code's mapping
        for code in icd10_codes:
            code = code.strip().strip("'\"")  # Remove quotes and whitespace
            if code not in icd10_to_traits:
                icd10_to_traits[code] = []
            icd10_to_traits[code].append(trait_idx)
    
    def get_trait_modules(individual_codes) -> List[int]:
        """Get trait modules for an individual's ICD10 codes."""
        if individual_codes is None or (isinstance(individual_codes, float) and pd.isna(individual_codes)):
            return []
        
        # Handle different input formats
        if isinstance(individual_codes, str):
            # Parse string representation if needed
            if individual_codes.startswith('['):
                codes = individual_codes.strip('[]').replace(' ', '').split(',')
                codes = [c.strip().strip("'\"") for c in codes if c.strip()]
            else:
                # Assume comma-separated or single code
                codes = [c.strip() for c in individual_codes.split(',')]
        elif isinstance(individual_codes, list):
            codes = individual_codes
        else:
            return []
        
        # Find all trait modules for this individual
        trait_modules: Set[int] = set()
        for code in codes:
            if code in icd10_to_traits:
                trait_modules.update(icd10_to_traits[code])
        
        return sorted(list(trait_modules))
    
    # Apply function to create trait modules for each individual
    result_df['trait_modules'] = result_df[icd10_column].swifter.apply(get_trait_modules)
    
    return result_df


# Example usage:
if __name__ == "__main__":
    # Example trait definitions (like your data)

    trait_df = pd.read_pickle("data/burden_test_trait_info.obj")
    
    # Example individual data
    individual_data = {
        'id': [1, 2, 3, 4, 5],
        'icd10_codes': [
            ['I48', 'I422'],  # Should match traits 22, 31, 61
            ['M069', 'G250'],  # Should match traits 91, 160
            ['I425'],  # Should match traits 31, 141
            ['A319', 'I999'],  # Should match trait 107, I999 not in any trait
            []  # No codes
        ]
    }
    individual_df = pd.DataFrame(individual_data)
    
    # Create trait modules
    result = create_trait_modules(trait_df.reset_index(), individual_df)
    print(result)
    print("\nTrait modules per individual:")
    for _, row in result.iterrows():
        print(f"Individual {row['id']}: {row['trait_modules']}")
