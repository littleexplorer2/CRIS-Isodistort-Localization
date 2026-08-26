"""Official Method 1 Conventional / Primitive lattice dropdown labels.

Source: saved ISODISTORT search page for EuAl4 Parent.cif (I4/mmm #139) in
``webpage_info/2. ISODISTORT_ search.html``. Labels and order match the
website ``isolattice`` / ``isoplattice`` selects.

Local enumeration yields the same lattice *classes* (unique under parent
point-group rotation + GL(3,Z)); these strings are the website's chosen
representatives and sort order for each class.
"""
from __future__ import annotations

# parent space-group number -> conventional / primitive option labels (official order)
METHOD1_LATTICE_OFFICIAL: dict[int, dict[str, list[str]]] = {
    139: {
        "conventional": [
            "(1,0,0),(0,1,0),(0,0,1)",
            "(1,1,0),(-1,1,0),(0,0,1)",
            "(1,-1,0),(1,1,0),(-1/2,1/2,1/2)",
            "(1,0,0),(0,1,0),(-1/2,-1/2,1/2)",
            "(1,0,-1),(0,-1,0),(-2,0,0)",
            "(2,0,0),(0,0,-2),(0,1,0)",
            "(2,2,0),(-2,2,0),(-1/2,-1/2,1/2)",
            "(-1/2,-1/2,1/2),(2,0,0),(0,2,0)",
            "(2,0,0),(0,2,0),(0,0,2)",
            "(2,-2,0),(2,2,0),(-1,1,1)",
            "(2,0,0),(0,2,0),(-1,-1,1)",
            "(-1,1,0),(-1,-1,0),(0,0,2)",
        ],
        "primitive": [
            "(-1/2,1/2,1/2),(1/2,-1/2,1/2),(1/2,1/2,-1/2)",
            "(-1/2,1/2,1/2),(1/2,-1/2,1/2),(1,1,0)",
            "(1,-1,0),(1,1,0),(0,0,1)",
            "(1,0,0),(0,1,0),(0,0,1)",
            "(1,0,1),(1/2,1/2,-1/2),(-1/2,1/2,1/2)",
            "(1,0,1),(1,0,-1),(0,1,0)",
            "(-1/2,-1/2,-3/2),(3/2,-1/2,1/2),(1/2,-3/2,-1/2)",
            "(-1,1,1),(1,-1,1),(1,1,-1)",
            "(0,-1,1),(0,1,1),(-1,0,-1)",
        ],
    },
}
