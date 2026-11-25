"""The portal: a read-only web view of a cgload database (D-111).

``cgload serve`` runs it; ``cgload index`` builds the search index it uses to
find features by name (D-112). The web parts need the ``portal`` extra
(``pip install "cgload[portal]"``); this package imports nothing from them at
import time, so ``cgload index`` works without them.
"""
