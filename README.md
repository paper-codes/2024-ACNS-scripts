# 2024-ACNS-scripts
Scripts used in the submission  [Quantum Circuit Design for the
    Lee-Brickell based Information Set Decoding](#authors-and-citations)

Source code for paper [Improving the Efficiency of Quantum Circuits for
Information Set Decoding](#authors-and-citations) submitted at . 

The bibtex can be downloaded [here](./bibtex.bib)

# Usage #
## Prange ##
```
# For PBP23
python -m measures.prange_launcher prange --out-format txt --expansion ours
# For Ess+21
python -m measures.prange_launcher esser21 --out-format txt
# For classic
python -m measures.leebrickell_launcher lee-brickell-classic --out-format txt -p 0
```


## Lee-Brickell ##

```
# For Ess+21
for p in {1..5}; do echo "#### $p ####" && python -m measures.leebrickell_launcher esser21-lb --out-format txt -p $p; done
# For PBP19 (Hybrid L-B)
for p in {1..5}; do echo "#### $p ####" && python -m measures.leebrickell_launcher lee-brickell-4 --out-format txt -p $p; done
# For classic
for p in {1..2}; do echo "#### $p ####" && python -m measures.leebrickell_launcher lee-brickell-classic --out-format txt -p $p; done
```


# Authors and citations #
The code here was used in the results of the following article

[PeBP23] _Perriello, Simone_ ; Barenghi, Alessandro ; Pelosi, Gerardo: /Improving the Efficiency of Quantum Circuits for
Information Set Decoding/. In: . — Citation Key:
[Accepted on December 2023]
[bibtex](/bibtex.bib)
