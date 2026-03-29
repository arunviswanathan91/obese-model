# Signature Curation with LLM

This script curates a custom immune gene signature database for PDAC and obesity research using Google Gemini as an AI assistant, with mandatory human review at every decision point.

**Script:** `Signature Curation with LLM.py`

The final signature database (2,143 signatures across 65 cell types, ~10 genes/signature) is browsable at [obese-pdac-model.streamlit.app](https://obese-pdac-model.streamlit.app) under the _Signature Explorer_ section.

---

## Overview

The curation pipeline runs in two sequential phases:

| Phase | Name          | Purpose                                                   |
| ----- | ------------- | --------------------------------------------------------- |
| 1     | Deduplication | Remove redundant signatures with high gene overlap        |
| 2     | Discovery     | Identify and add biologically relevant missing signatures |

---

## Phase 1 — Deduplication

All signature pairs within each cell type were compared by **Overlap Coefficient**:

> O = |A ∩ B| / min(|A|, |B|)

Pairs with overlap > 0.50 were flagged and reviewed. For each flagged pair, Gemini provided a recommendation (keep one, keep both, or merge), and the researcher made the final call from five options: keep signature 1, keep signature 2, keep both, merge into signature 1, or merge into signature 2.

> **Note:** only the Overlap Coefficient is computed in the code.

---

## Phase 2 — Discovery

For each cell type, Gemini was prompted to suggest up to **3 missing signatures** relevant to obesity-induced dysfunction in PDAC. Each suggestion was validated programmatically before being shown for human review.

**Acceptance criteria for new signatures:**

| Criterion                                   | Value                         |
| ------------------------------------------- | ----------------------------- |
| Gene count                                  | 8–12 genes                    |
| Maximum overlap with any existing signature | < 20% (Overlap Coefficient)   |
| Human approval                              | Required (y/n per suggestion) |

---

## LLM Configuration

| Setting            | Value                                                                                  |
| ------------------ | -------------------------------------------------------------------------------------- |
| API                | Google Gemini (google-generativeai)                                                    |
| Model selection    | Auto — priority: gemini-1.5-flash → gemini-2.0-flash-exp → gemini-1.5-pro → gemini-pro |
| Context            | Pancreatic Cancer (PDAC) + Obesity                                                     |
| Response format    | JSON only                                                                              |
| Retries on failure | 3 (with backoff)                                                                       |

---

## Human-in-the-Loop

Every AI recommendation — both deduplication decisions and new signature additions — required explicit researcher approval before being applied. All actions were written to a session log (`Interactive_Log.txt`).

### Interactive Prompts

**Phase 1 — Deduplication** (one prompt per flagged pair):

The script displays both signature names, the overlap percentage, and a gene breakdown (shared genes, unique to each), then shows the AI recommendation and its reasoning before asking:

```text
Your Decision (1-5, S):
  1 — Keep Signature 1 only (delete Signature 2)
  2 — Keep Signature 2 only (delete Signature 1)
  3 — Keep both
  4 — Merge into Signature 1
  5 — Merge into Signature 2
  S — Skip
```

**Phase 2 — Discovery** (one prompt per valid suggestion):

The script displays the proposed signature name, biological rationale, and gene list before asking:

```text
Add this signature? [y/n]:
```
