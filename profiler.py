import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

# Step 2: Fetch data from NCBI
Entrez.email = "your_email@example.com"  # Replace with your real email


def fetch_gene_record(accession_id: str):
  print(f"[1/4] Fetching {accession_id} from NCBI...")
  with Entrez.efetch(
      db="nucleotide", id=accession_id, rettype="gb", retmode="text"
  ) as handle:
    record = SeqIO.read(handle, "genbank")
  return record


# Step 3: Extract CDS and Translate
def get_translated_protein(seq_record):
  print("[2/4] Translating coding sequence to protein...")
  for feature in seq_record.features:
    if feature.type == "CDS":
      cds_seq = feature.extract(seq_record.seq)
      protein_seq = cds_seq.translate(to_stop=True)
      return cds_seq, protein_seq
  return seq_record.seq, seq_record.seq.translate(to_stop=True)


# Step 4: Scan for Stress Motifs
def scan_motifs(sequence: Seq, motifs: list[str]) -> dict[str, list[int]]:
  print("[3/4] Scanning for stress-responsive motifs...")
  seq_str = str(sequence).upper()
  motif_locations = {}
  for motif in motifs:
    positions = []
    start = 0
    while True:
      pos = seq_str.find(motif, start)
      if pos == -1:
        break
      positions.append(pos)
      start = pos + 1
    motif_locations[motif] = positions
  return motif_locations


# Step 5: Plot GC Content
def plot_gc_sliding_window(
    sequence: Seq,
    gene_name: str,
    window_size: int = 50,
    step: int = 10,
    save_path: str = None,
):
  print("[4/4] Generating GC-content profile plot...")
  seq_str = str(sequence).upper()
  gc_values = []
  positions = []

  for i in range(0, len(seq_str) - window_size + 1, step):
    subseq = seq_str[i : i + window_size]
    gc = (subseq.count("G") + subseq.count("C")) / window_size * 100
    gc_values.append(gc)
    positions.append(i + (window_size // 2))

  plt.figure(figsize=(10, 4))
  plt.plot(positions, gc_values, color="forestgreen", lw=1.5)
  avg_gc = sum(gc_values) / len(gc_values)
  plt.axhline(
      avg_gc,
      color="red",
      linestyle="--",
      label=f"Avg GC ({avg_gc:.1f}%)",
  )
  plt.title(f"GC-Content Profile - {gene_name}")
  plt.xlabel("Nucleotide Position (bp)")
  plt.ylabel("GC Content (%)")
  plt.grid(True, linestyle=":", alpha=0.6)
  plt.legend()
  plt.tight_layout()

  if save_path:
    plt.savefig(save_path, dpi=300)
    print(f"Plot saved successfully to {save_path}")
  plt.show()


# Execution Flow
if __name__ == "__main__":
  # 1. Download
  gene_id = "NM_118680"  # Arabidopsis DREB1A
  record = fetch_gene_record(gene_id)
  print(f"Retrieved: {record.description[:60]}... ({len(record.seq)} bp)\n")

  # 2. Translate
  cds, protein = get_translated_protein(record)
  print(f"CDS Length: {len(cds)} bp | Protein Length: {len(protein)} aa")
  print(f"First 25 amino acids: {protein[:25]}\n")

  # 3. Motif Detection
  motifs = ["GCCGAC", "ACCGAC"]
  matches = scan_motifs(record.seq, motifs)
  for m, hits in matches.items():
    print(f"Motif '{m}': {len(hits)} occurrence(s) at positions {hits}")
  print()

  # 4. Visualization
  plot_gc_sliding_window(
      record.seq, gene_name="Arabidopsis DREB1A", save_path="dreb1a_gc_plot.png"
  )