from molop.io.base_models.SearchPattern import MolOPPattern


class FCHKPatterns:
    NMR_ROUTE = MolOPPattern(
        content_pattern=r"(?i)\bNMR\s*=\s*(?P<options>\([^)]*\)|[^\s]+)",
        description="Gaussian NMR route options stored in a formatted checkpoint.",
    )


fchk_patterns = FCHKPatterns()
