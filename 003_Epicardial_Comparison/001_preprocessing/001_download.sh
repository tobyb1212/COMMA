#!/usr/bin/env bash
cd /COMMA/results/scrna/downstream_analysis/epicardial_coronary/data/
curl --location --fail https://service.azul.data.humancellatlas.org/manifest/files/ksQwlKVkY3A1OaRjdXJsxBDaz5YIyS5TWoGEjOWzk1v8xBD96PWmHcdWv4QaD817EDL0xCB89-wGFPlsECvscQYnltFw5P0FHi8pQdSCKkeSAYAqdg | curl --retry 15 --retry-delay 10 --config -
mv */* .