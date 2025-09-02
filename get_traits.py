import pandas as pd

def main():
    # read .obj file
    df = pd.read_pickle('data/burden_test_trait_info.obj')
    print(df.head())
    print(df.columns)
    print(df.icd10)

if __name__ == "__main__":
    main()