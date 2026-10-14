# Slice 3 break-once proof (2026-10-15)

Registration-slice-3.md section 6, gate 3. `tools/slice-3/mutants.py` (slice 1's runner) on a clone of `9033206`, one
mutant at a time, `swift test -c release --filter 'CTFTests|CTFTruthTests'`.

| Mutant | Rule broken | Killed by |
|---|---|---|
| C01 | E_t without the square | envelopesEqualTheirClosedForms, hrtemImagesEqualAbTEM, kernelsEqualAbTEM |
| C02 | E_t with λ multiplied | envelopesEqualTheirClosedForms, hrtemImagesEqualAbTEM, kernelsEqualAbTEM |
| C03 | β not converted from mrad | envelopesEqualTheirClosedForms, hrtemImagesEqualAbTEM, kernelsEqualAbTEM |
| C04 | ∂χ/∂φ's first ½ dropped | hrtemImagesEqualAbTEM, kernelsEqualAbTEM |
| C05 | ∂χ/∂α missing C30 | envelopesEqualTheirClosedForms, hrtemImagesEqualAbTEM, kernelsEqualAbTEM |
| C06 | sign(β) dropped (equivalent for β ≥ 0) | equivalent: the port refuses β < 0, so sign(β) is 1 or 0 and dropping it changes nothing |
| C07 | the soft aperture made hard | hrtemImagesEqualAbTEM, kernelsEqualAbTEM |
| C08 | the kernel conjugated | hrtemImagesEqualAbTEM, kernelsEqualAbTEM, weakPhaseContrastAtScherzer |
| C09 | the temporal envelope applied twice | hrtemImagesEqualAbTEM, kernelsEqualAbTEM |
| C10 | apply_ctf scaled wrongly (×2 on the real part) | hrtemImagesEqualAbTEM, identityCTFLeavesTheImage, vacuumImageIsFlat, weakPhaseContrastAtScherzer |
| C11 | Scherzer's 3/2 replaced by 4/3 | closedFormsEqualAbTEMAndTheTextbook |
| C12 | point resolution with a cube root | closedFormsEqualAbTEMAndTheTextbook |
| C13 | Nyquist with 2 instead of 4 | closedFormsEqualAbTEMAndTheTextbook |
| C14 | a negative angular spread accepted | refusals |
| C15 | the aperture not applied | hrtemImagesEqualAbTEM, kernelsEqualAbTEM, weakPhaseContrastAtScherzer |

Result: 14 of 15 killed, C06 equivalent (planned as such in section 6). Every slice-3 test function appears in the last column.
