"""
Quick demo: verifies the cBioPortal connection and pulls key data from the
Pediatric Wilms' Tumor (TARGET, 2018) dataset.

Run with:
    python demo.py
"""

from wilms_tumor import WilmsTumorDataset

def main():
    wt = WilmsTumorDataset()

    # ── 1. Study summary ─────────────────────────────────────────────
    wt.summary()

    # ── 2. Clinical data ─────────────────────────────────────────────
    print("\n[1] Fetching patient clinical data...")
    clinical = wt.clinical_data()
    print(f"    Shape: {clinical.shape}")
    print(f"    Columns: {list(clinical.columns)}")
    if "OS_MONTHS" in clinical.columns:
        print(f"    OS_MONTHS stats:\n{clinical['OS_MONTHS'].astype(float, errors='ignore').describe()}")

    # ── 3. Top mutated genes ─────────────────────────────────────────
    print("\n[2] Top 15 most frequently mutated genes (all 657 samples)...")
    freq = wt.mutation_frequency(top_n=15)
    print(freq.to_string())

    # ── 4. Mutations for key Wilms' genes ────────────────────────────
    # WTX was renamed to AMER1 in HGNC — use AMER1
    wilms_genes = ["WT1", "CTNNB1", "TP53", "AMER1", "SIX1", "SIX2", "DROSHA", "DGCR8"]
    print(f"\n[3] Mutations for key Wilms' Tumor genes: {wilms_genes}")
    for gene in wilms_genes:
        try:
            df = wt.mutations_for_gene(gene)
            n_samples = df["sampleId"].nunique() if not df.empty else 0
            print(f"    {gene:<10}: {n_samples:>3} mutated samples")
        except Exception as e:
            print(f"    {gene:<10}: error — {e}")

    # ── 5. CNA overview ──────────────────────────────────────────────
    print("\n[4] Fetching CNA data (GISTIC, 129 samples)...")
    cna = wt.cna()
    print(f"    Shape: {cna.shape}")
    if not cna.empty:
        val_col = "value" if "value" in cna.columns else cna.columns[-1]
        print(f"    CNA value distribution:\n{cna[val_col].value_counts().sort_index()}")

    # ── 5. Expression for canonical Wilms' genes ─────────────────────
    print("\n[5] Fetching RNA-Seq expression for canonical genes...")
    expr = wt.expression_for_genes(["WT1", "CTNNB1", "TP53", "IGF2"])
    print(f"    Shape: {expr.shape}")
    if not expr.empty:
        print(f"    Preview:\n{expr.head(5).to_string(index=False)}")

    print("\nDone. Connection verified.")

if __name__ == "__main__":
    main()
