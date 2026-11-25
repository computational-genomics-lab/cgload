#!/usr/bin/env python3
"""Figure 3: three pages of the cgload portal, photographed in a browser.

    python3 paper/screenshots.py [--out paper]

The genome is the repository's own test record of the Dscam1 region of
*Drosophila melanogaster* (tests/fixtures/genbank_dscam1_CG12164.gbff, NCBI
NT_033778.4), loaded with the record's own sequence, so every base and every
protein the pages show is real. The script loads it into a new SQLite
database, indexes it, serves it with `cgload serve`, opens three pages in
Chromium (Playwright) and puts them side by side in OUT/figures/fig3.pdf and
.png. Needs the [portal] extra, playwright with its Chromium, and matplotlib.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import socket
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.image import imread  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
FIXTURE = REPO / "tests" / "fixtures" / "genbank_dscam1_CG12164.gbff"
FULL_WIDTH = 170 / 25.4      # inches: the journal's full page width
VIEW = 1280                  # CSS pixels across each page
SCALE = 2                    # device pixels per CSS pixel: 2,560 across 166 mm, 392 dpi
INK = "#0b0b0b"


def run(*args: str, env: dict) -> None:
    # This Python and cgload's own commands, with paths made here.
    subprocess.run([sys.executable, "-m", "cgload", *args], check=True, env=env,  # noqa: S603
                   stdout=subprocess.DEVNULL)


def write_fasta(gbff: Path, target: Path) -> None:
    """The record's own sequence, from its ORIGIN section, as FASTA."""
    text = gbff.read_text()
    seqid = re.search(r"^VERSION\s+(\S+)", text, re.M).group(1)
    bases = "".join(re.sub(r"[^a-z]", "", line)
                    for line in text.split("\nORIGIN", 1)[1].splitlines()).upper()
    with target.open("w") as out:
        out.write(f">{seqid}\n")
        for start in range(0, len(bases), 60):
            out.write(bases[start:start + 60] + "\n")


def build(folder: Path, env: dict) -> tuple[str, Path]:
    url = f"sqlite:///{folder / 'portal.db'}"
    data = folder / "data"
    fly = folder / "fly"
    fly.mkdir()
    (fly / FIXTURE.name).write_text(FIXTURE.read_text())
    write_fasta(FIXTURE, fly / "region.fa")
    (fly / "organism.toml").write_text(
        '[organism]\ngenus = "Drosophila"\nspecies = "melanogaster"\n'
        'strain = "reference"\nncbi_taxon_id = 7227\n\n'
        '[assembly]\nname = "Release 6, Dscam1 region"\naccession = "NT_033778.4"\n'
        'fasta = "region.fa"\n\n'
        f'[[files.annotation]]\npath = "{FIXTURE.name}"\n')
    run("init", "--db-url", url, env=env)
    run("load", "--config", str(fly / "organism.toml"), "--db-url", url, "--data-dir",
        str(data), env=env)
    run("index", "--db-url", url, env=env)
    return url, data


def free_port() -> int:
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        return probe.getsockname()[1]


def get_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=10) as response:  # noqa: S310 - localhost
        return json.loads(response.read())


def photograph(base: str, folder: Path) -> list[Path]:
    from playwright.sync_api import sync_playwright

    hits = get_json(f"{base}/api/v1/genomes/1/search?q=Dscam1")["results"]
    gene = next(hit["feature"]["row"] for hit in hits if hit["feature"]["type"] == "gene")
    shots = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport={"width": VIEW, "height": 900},
                                device_scale_factor=SCALE)
        outside = []
        page.on("request", lambda r: outside.append(r.url) if not r.url.startswith(base)
                else None)

        def shoot(path: str, top: str, bottom: str, name: str, *, wait_igv=False) -> None:
            page.goto(base + path, wait_until="networkidle")
            if wait_igv:
                page.wait_for_selector("#igv canvas", timeout=20000)
                page.wait_for_timeout(2500)
            first = page.locator(top).first.bounding_box()
            last = page.locator(bottom).first.bounding_box()
            y = max(first["y"] - 12, 0)
            clip = {"x": 0, "y": y, "width": VIEW,
                    "height": last["y"] + last["height"] - y + 2}
            target = folder / f"{name}.png"
            page.screenshot(path=str(target), clip=clip, full_page=True)
            if wait_igv:
                box = page.locator("#igv").bounding_box()
                close_up(target, (box["x"] - clip["x"], box["y"] - y, box["width"],
                                  box["height"]))
                # The key under the gene model; the figure's legend says it instead.
                key = page.locator("text=Boxes are exons").first.bounding_box()
                leave_out(target, key["y"] - y - 6, key["y"] + key["height"] - y + 6)
            shots.append(target)

        # a: a gene, its transcripts drawn to scale and the genome browser at its locus
        shoot(f"/g/1/f/{gene}", "header", "#igv", "a", wait_igv=True)
        # b: a search by words of a product, with how each feature matched
        shoot("/g/1/search?q=cell+adhesion+molecule", "main p:has-text('features match')",
              ".pager", "b")
        # c: the genome's page: how it was loaded, checks and input files
        shoot("/g/1/", "h2:has-text('How this genome was loaded')",
              "h2:has-text('How this genome was loaded') ~ table >> nth=1", "c")
        browser.close()
    if outside:
        raise SystemExit(f"a page asked another site for {outside[:3]}")
    return shots


def close_up(image_path: Path, box: tuple[float, float, float, float]) -> None:
    """Take out the empty rows at the foot of the genome browser's track area, so
    the figure shows the tracks without the room the page leaves below them."""
    import numpy as np

    image = imread(str(image_path))
    x, y, w, h = (round(value * SCALE) for value in box)
    # The tracks, without the browser's border, its scroll bar and its gear icons.
    area = image[y:y + h - round(4 * SCALE), x + round(12 * SCALE):x + w - round(70 * SCALE), :3]
    inked = np.nonzero((area < 0.85).any(axis=2).any(axis=1))[0]
    last = y + int(inked.max()) if inked.size else y + h
    cut_from, cut_to = last + round(14 * SCALE), y + h - round(3 * SCALE)
    if cut_to > cut_from:
        image = np.concatenate([image[:cut_from], image[cut_to:]])
        plt.imsave(str(image_path), image)


def leave_out(image_path: Path, top: float, bottom: float) -> None:
    """Take the rows between ``top`` and ``bottom`` (CSS pixels) out of an image."""
    import numpy as np

    image = imread(str(image_path))
    image = np.concatenate([image[:round(top * SCALE)], image[round(bottom * SCALE):]])
    plt.imsave(str(image_path), image)


def compose(shots: list[Path], out: Path) -> None:
    images = [imread(str(path)) for path in shots]
    gap, label = 0.12, 0.16                     # inches
    width = FULL_WIDTH - label                  # the pages; the letters take the rest
    heights = [width * image.shape[0] / image.shape[1] for image in images]
    total = sum(heights) + gap * (len(images) - 1)
    figure = plt.figure(figsize=(FULL_WIDTH, total))
    y = total
    for letter, image, height in zip("abc", images, heights, strict=True):
        y -= height
        axis = figure.add_axes((label / FULL_WIDTH, y / total, width / FULL_WIDTH,
                                height / total))
        axis.imshow(image, interpolation="lanczos")
        axis.set_axis_off()
        figure.text(0, (y + height) / total, letter, ha="left", va="top", fontsize=9,
                    fontweight="bold", color=INK)
        y -= gap
    out.mkdir(parents=True, exist_ok=True)
    figure.savefig(out / "fig3.pdf", dpi=375, metadata={"CreationDate": None})
    figure.savefig(out / "fig3.png", dpi=300, metadata={"Software": None})
    plt.close(figure)
    print(f"Figure 3: {FULL_WIDTH * 25.4:.0f} x {total * 25.4:.0f} mm")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=HERE)
    args = parser.parse_args(argv)
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env.get("PYTHONPATH", "")
    env.pop("CGLOAD_DB_URL", None)
    with tempfile.TemporaryDirectory() as scratch:
        folder = Path(scratch)
        url, data = build(folder, env)
        port = free_port()
        server = subprocess.Popen(  # noqa: S603 - this Python and cgload serve
            [sys.executable, "-m", "cgload", "serve", "--db-url", url, "--data-dir",
             str(data), "--port", str(port), "--title", "cgload portal"],
            env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        base = f"http://127.0.0.1:{port}"
        try:
            for _ in range(100):
                try:
                    get_json(base + "/healthz")
                    break
                except OSError:
                    time.sleep(0.2)
            shots = photograph(base, folder)
            compose(shots, args.out / "figures")
        finally:
            server.terminate()
            server.wait(timeout=30)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
