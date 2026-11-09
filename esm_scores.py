"""Opt-in protein-language-model scoring (ESM-1v / ESM-2 zero-shot).

Meier et al., NeurIPS 2021 showed that a masked-LM log-likelihood ratio is
state-of-the-art zero-shot for mutation effects -- comparable to supervised
tools and complementary to FoldX, with no local binary beyond the model
weights.  This adapter is intentionally opt-in: it activates only when the
`ECD_ESM_MODEL` environment variable is set AND torch + fair-esm import
cleanly.  Anything else degrades to None and the run states that the PLM
evidence line was unavailable; no other behaviour changes.
"""
from __future__ import annotations

import math
import os
from typing import Any

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
_ACTIVATION_ENV = "ECD_ESM_MODEL"


class PlmScorer:
    """Lazy ESM masked-marginal scorer; construct via :func:`open_scorer`."""

    def __init__(self, model_name: str, model, alphabet, batch_converter) -> None:
        self.model_name = model_name
        self._model = model
        self._alphabet = alphabet
        self._batch_converter = batch_converter
        self._model.eval()

    def substitution_llrs(self, sequence: str) -> dict[tuple[int, str], float]:
        """LLR of every substitution at every position, in natural log units."""
        import torch

        data = [("protein", sequence)]
        _, _, tokens = self._batch_converter([data[0]])
        with torch.no_grad():
            logits = self._model(tokens)["logits"][0]  # (L+2, vocab)
        for_spec = self._alphabet.get_idx("<mask>")
        results: dict[tuple[int, str], float] = {}
        for index, wt in enumerate(sequence, start=1):
            masked = tokens[0].clone()
            masked[index] = for_spec
            with torch.no_grad():
                masked_logits = self._model(masked.unsqueeze(0))["logits"][0, index]
            log_probs = torch.log_softmax(masked_logits, dim=-1)
            wt_idx = self._alphabet.get_idx(wt)
            wt_ll = float(log_probs[wt_idx])
            for mutant in AMINO_ACIDS:
                if mutant == wt:
                    continue
                ll = float(log_probs[self._alphabet.get_idx(mutant)])
                results[(index, mutant)] = round(ll - wt_ll, 4)
        del logits
        return results


def open_scorer() -> tuple[PlmScorer, None] | tuple[None, str]:
    """Return (scorer, None) or (None, reason); never raises."""
    model_name = os.environ.get(_ACTIVATION_ENV, "").strip()
    if not model_name:
        return None, f"{_ACTIVATION_ENV} is not set; PLM scoring stays disabled"
    try:
        import esm
        import torch  # noqa: F401
    except Exception as exc:
        return None, f"torch/fair-esm unavailable ({exc.__class__.__name__}); PLM scoring disabled"
    try:
        model, alphabet = getattr(esm, model_name)()
        batch_converter = alphabet.get_batch_converter()
        return PlmScorer(model_name, model, alphabet, batch_converter), None
    except Exception as exc:
        return None, f"could not load {model_name} ({exc}); PLM scoring disabled"


def score_choices(choices: list[Any], sequence: str, ecd_start: int) -> dict[tuple[int, str], float] | None:
    """Convenience wrapper; returns None (with no work done) when unusable."""
    scorer, _reason = open_scorer()
    if scorer is None:
        return None
    if len(sequence) > 1200:
        return None
    llrs = scorer.substitution_llrs(sequence)
    return {
        (position + ecd_start - 1, mutant): value
        for (position, mutant), value in llrs.items()
    }


def llr_to_flag(llr: float | None) -> str:
    """Consensus flag for a choice's reasons list; None -> unavailable note."""
    if llr is None:
        return "PLM evidence unavailable"
    if llr >= -0.5:
        return "PLM tolerated"
    if llr >= -1.5:
        return "PLM mildly deleterious"
    if llr >= -3.0:
        return "PLM deleterious"
    return "PLM strongly deleterious"


def plausible_llr(value: float) -> bool:
    return isinstance(value, float) and math.isfinite(value)
