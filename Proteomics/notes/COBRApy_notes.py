from __future__ import annotations
from collections import defaultdict
from typing import Dict, Set, Iterable, Optional

import cobra


def _norm_subsystem(s: str) -> str:
    s = (s or "").strip()
    return "Unassigned" if s == "" else s


def _get_reaction_subsystems(rxn: cobra.Reaction) -> Set[str]:
    """
    Robustly extract subsystem(s) from a COBRApy Reaction.
    Human1 commonly uses rxn.subsystem (string), but some models store lists or notes/annotation.
    """
    subs: Set[str] = set()

    # 1) Most common: rxn.subsystem (string)
    if hasattr(rxn, "subsystem"):
        val = getattr(rxn, "subsystem")
        if isinstance(val, str) and val.strip():
            subs.add(_norm_subsystem(val))
        elif isinstance(val, (list, tuple, set)):
            for x in val:
                if isinstance(x, str) and x.strip():
                    subs.add(_norm_subsystem(x))

    # 2) Sometimes: rxn.subsystems (list-like)
    if hasattr(rxn, "subsystems"):
        val = getattr(rxn, "subsystems")
        if isinstance(val, str) and val.strip():
            subs.add(_norm_subsystem(val))
        elif isinstance(val, (list, tuple, set)):
            for x in val:
                if isinstance(x, str) and x.strip():
                    subs.add(_norm_subsystem(x))

    # 3) Some reconstructions: rxn.notes may contain 'subsystem'
    if hasattr(rxn, "notes") and isinstance(rxn.notes, dict):
        for k in ("subsystem", "Subsystem", "SUBSYSTEM"):
            if k in rxn.notes:
                val = rxn.notes[k]
                if isinstance(val, str) and val.strip():
                    subs.add(_norm_subsystem(val))
                elif isinstance(val, (list, tuple, set)):
                    for x in val:
                        if isinstance(x, str) and x.strip():
                            subs.add(_norm_subsystem(x))

    # 4) Or rxn.annotation
    if hasattr(rxn, "annotation") and isinstance(rxn.annotation, dict):
        for k in ("subsystem", "Subsystem", "pathway", "Pathway"):
            if k in rxn.annotation:
                val = rxn.annotation[k]
                if isinstance(val, str) and val.strip():
                    subs.add(_norm_subsystem(val))
                elif isinstance(val, (list, tuple, set)):
                    for x in val:
                        if isinstance(x, str) and x.strip():
                            subs.add(_norm_subsystem(x))

    if not subs:
        subs.add("Unassigned")
    return subs


def build_metabolite_subsystem_maps(
    model: cobra.Model,
    *,
    exclude_currency_mets: bool = True,
    currency_met_ids: Optional[Iterable[str]] = None,
    min_mets_per_subsystem: int = 5,
    drop_unassigned: bool = True,
) -> tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    """
    Map metabolites to subsystems via participating reactions.
    Returns:
      met_to_subs: met.id -> set(subsystem names)
      sub_to_mets: subsystem name -> set(met.id)
    """
    if currency_met_ids is None:
        # You can expand this list; keep it conservative by default.
        currency_met_ids = {
            "h2o", "h", "pi", "ppi", "co2", "o2", "nh4", "na1", "k",
            "atp", "adp", "amp", "nad", "nadh", "nadp", "nadph",
            "coa", "accoa", "hco3", "so4", "cl"
        }

    currency_met_ids = set(currency_met_ids)

    met_to_subs: Dict[str, Set[str]] = defaultdict(set)
    sub_to_mets: Dict[str, Set[str]] = defaultdict(set)

    # COBRApy contract: model.reactions iterable, rxn.metabolites dict of met->coeff
    for rxn in model.reactions:
        subs = _get_reaction_subsystems(rxn)
        for met in rxn.metabolites.keys():
            mid = met.id
            if exclude_currency_mets and (mid in currency_met_ids):
                continue
            for s in subs:
                if drop_unassigned and s == "Unassigned":
                    continue
                met_to_subs[mid].add(s)
                sub_to_mets[s].add(mid)

    # Filter tiny subsystems (helps PSI stability)
    if min_mets_per_subsystem > 1:
        keep_subs = {s for s, mets in sub_to_mets.items() if len(mets) >= min_mets_per_subsystem}
        sub_to_mets = {s: mets for s, mets in sub_to_mets.items() if s in keep_subs}
        met_to_subs = {
            mid: {s for s in subs if s in keep_subs}
            for mid, subs in met_to_subs.items()
            if any(s in keep_subs for s in subs)
        }

    return met_to_subs, sub_to_mets


if __name__ == "__main__":
    # Example usage
    model = cobra.io.load_matlab_model("Human1.mat")  # or read_sbml_model(...)
    met_to_subs, sub_to_mets = build_metabolite_subsystem_maps(model)

    print("Example metabolite -> subsystems:")
    for mid in list(met_to_subs.keys())[:5]:
        print(mid, sorted(list(met_to_subs[mid]))[:5])

    print("\nTop 5 largest subsystems:")
    top = sorted(sub_to_mets.items(), key=lambda kv: len(kv[1]), reverse=True)[:5]
    for s, mets in top:
        print(s, len(mets))
