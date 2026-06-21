# Data

## old_faithful_waiting.csv
Classic Old Faithful Geyser waiting-time dataset (Azzalini & Bowman, 1990).
272 consecutive eruption waiting times (minutes) from August 1985.
Public domain; available in R (`faithful$waiting`) and many statistics textbooks.

**Columns:** `observation` (sequential index), `waiting_time_min` (minutes between eruptions)

## dsdp609_hsg_record.csv
Hematite-Stained Grain (HSG) percentage record from DSDP Site 94-609,
North Atlantic (49.88°N, 24.24°W). Original counts by Bond et al. (1999),
re-examined by Obrochta et al. (2012). 310 measurements spanning MIS 4–2
(14.4–70.8 ka BP).

**Source:** PANGAEA, https://doi.org/10.1594/PANGAEA.834688
**License:** CC-BY-3.0

**Columns:**
- `depth_mbsf`: depth in meters below seafloor
- `age_ka_bp`: age in thousands of years before present (2000 CE)
- `hsg_pct`: hematite-stained grains (% of lithics, 63–150 µm)
- `volcanic_glass_pct`: Icelandic glass (% of lithics)
- `detrital_carbonate_pct`: detrital carbonate (% of lithics)

## ird_peaks_extracted.csv
35 HSG peaks extracted algorithmically from the DSDP 609 record using
`scipy.signal.find_peaks` with thresholds: HSG ≥ 15%, minimum inter-peak
distance = 3 samples, prominence ≥ 4 percentage points. These peaks mark
ice-rafted debris pulses (Bond events and Heinrich events).

**Columns:**
- `peak_number`: sequential peak index
- `age_ka_bp`: age of peak in ka BP
- `hsg_pct`: HSG percentage at peak
- `spacing_to_next_ka`: interval to next peak in kiloyears (blank for last peak)
