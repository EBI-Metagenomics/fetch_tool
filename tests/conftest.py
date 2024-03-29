#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2018-2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import os


def pytest_generate_tests(metafunc):
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    os.environ["FETCH_TOOL_CONFIG"] = os.path.join(root, "config/testing.json")
    # Follow the instructions on the README file to install it
    os.environ.setdefault("ASPERA_BIN", f"{root}/aspera-cli/cli/bin/ascp")
    os.environ.setdefault(
        "ASPERA_CERT", f"{root}/aspera-cli/cli/etc/asperaweb_id_dsa.openssh"
    )
