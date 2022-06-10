./misc/makeCondorJobForAnalysis.py  \
    bin/b0ToDstarSMinus.exe \
    configs/bsToDsDminusMuPMuM.tpl.cfg \
    results/BsToDstarXX_CMS2 \
    10 \
    50000 \
    bs2DStarCMS2
sleep 2

./misc/makeCondorJobForAnalysis.py  \
    bin/bsToDstarSMinusFullAcceptance.exe \
    configs/bsToDsDminusMuPMuM.tpl.cfg \
    results/BsToDstarXX_FA2 \
    10 \
    50000 \
    bs2DStarFA2
sleep 2

./misc/makeCondorJobForAnalysis.py  \
    bin/b0ToDstarSMinus.exe \
    configs/b0ToDDminusMuPMuM.tpl.cfg \
    results/B0ToDstarXX_CMS2 \
    10 \
    50000 \
    b02DStarCMS2
sleep 2

./misc/makeCondorJobForAnalysis.py  \
    bin/b0ToDstarSMinusFullAcceptance.exe \
    configs/b0ToDDminusMuPMuM.tpl.cfg \
    results/B0ToDstarXX_FA2 \
    10 \
    50000 \
    b02DStarFA2
sleep 2
