/* cgload portal: starts the genome browser where a page has one.
   Loaded as a file of the portal itself: the pages carry no inline script. */
(function () {
  "use strict";
  function start() {
    var element = document.getElementById("igv");
    if (!element || typeof igv === "undefined") {
      return;
    }
    var data = element.dataset;
    var options = {
      reference: { id: data.name, name: data.name, fastaURL: data.fasta, indexURL: data.fai },
      locus: data.locus,
      loadDefaultGenomes: false,
      showSampleNames: false,
      search: { url: data.locate, type: "json", coords: 1, timeout: 5000 },
      tracks: []
    };
    if (data.gff && data.tbi) {
      options.tracks.push({
        name: "Annotation, as stored by cgload",
        type: "annotation",
        format: "gff3",
        url: data.gff,
        indexURL: data.tbi,
        displayMode: "EXPANDED",
        height: 320,
        visibilityWindow: 2000000
      });
    }
    igv.createBrowser(element, options).catch(function (error) {
      element.textContent = "The genome browser could not start: " + error;
    });
  }
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start);
  } else {
    start();
  }
})();
