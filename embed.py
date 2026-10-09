"""DOC-2-072: mean-pooled ESM-2 35M embeddings -> emb.npy (row order = dataset.tsv)."""
import numpy as np, pandas as pd, torch, time
from transformers import AutoTokenizer, EsmModel
torch.set_num_threads(2); d = pd.read_csv("dataset.tsv", sep="\t")
tok = AutoTokenizer.from_pretrained("facebook/esm2_t12_35M_UR50D"); m = EsmModel.from_pretrained("facebook/esm2_t12_35M_UR50D").eval()
order = np.argsort(d.length.values, kind="stable"); emb = np.zeros((len(d), m.config.hidden_size), dtype=np.float32); t0 = time.time(); B = 16
with torch.no_grad():
    for i in range(0, len(d), B):
        idx = order[i:i + B]; b = tok([d.sequence.iloc[j] for j in idx], return_tensors="pt", padding=True)
        h = m(**b).last_hidden_state; mask = b["attention_mask"].clone(); mask[:, 0] = 0
        for k, j in enumerate(idx):
            L = int(b["attention_mask"][k].sum()); mask[k, L - 1] = 0
        w = mask.unsqueeze(-1).float(); emb[idx] = ((h * w).sum(1) / w.sum(1)).numpy()
        if (i // B) % 20 == 0: print(i, len(d), round(time.time() - t0), flush=True)
np.save("emb.npy", emb); print("EMBED_DONE", emb.shape)
