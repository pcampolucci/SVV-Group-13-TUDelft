#!/bin/bash

printf "***************************************************"
printf "\e[5m\n           RUNNING TESTS           \n\e[25m"
printf "***************************************************"

printf "\n\nRunning tests on functions calculating geometry\n\n"
pytest geometry_tests.py

printf "\n\nRunning tests on functions calculating loadings\n\n"
pytest loading_tests.py

printf "\n\nRunning tests on infrastructure robustness\n\n"
pytest infra_tests.py