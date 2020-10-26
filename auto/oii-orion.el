(TeX-add-style-hook
 "oii-orion"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("mnras" "useAMS" "usenatbib")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("babel" "spanish" "es-minimal" "english") ("inputenc" "utf8") ("newtxmath" "stix2" "smallerops") ("enumitem" "shortlabels")))
   (TeX-run-style-hooks
    "latex2e"
    "mnras"
    "mnras10"
    "babel"
    "inputenc"
    "graphicx"
    "xcolor"
    "hyperref"
    "siunitx"
    "newtxtext"
    "newtxmath"
    "booktabs"
    "etoolbox"
    "enumitem")
   (TeX-add-symbols
    '("chem" 1)
    '("Wav" 1)
    '("ION" 2)
    "hii"
    "nii"
    "oiii"
    "oii"
    "Fion"
    "ionpar")
   (LaTeX-add-labels
    "firstpage"
    "sec:introduction"
    "sec:conclusions"
    "lastpage")
   (LaTeX-add-bibliographies
    "BibdeskLibrary")
   (LaTeX-add-counters
    "ionstage"))
 :latex)

