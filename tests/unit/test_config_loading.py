#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2018-2024 EMBL - European Bioinformatics Institute
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

from fetchtool import fetch_reads

FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "fixtures"))


class TestConfigLoading:
    def test_config_defaults(self):
        fetch = fetch_reads.FetchReads(argv=["-p", "ERP001736"])
        assert fetch.config == {
            "ena_api_username": "",
            "ena_api_password": "",
            "url_max_attempts": 5,
            "fire_endpoint": "https://hl.fire.sdo.ebi.ac.uk",
            "fire_ena_bucket": "era-public",
            "fire_access_key_id": "",
            "fire_secret_access_key": "",
        }

    def test_config_override_with_json_file(self):
        fetch = fetch_reads.FetchReads(argv=["-p", "ERP001736", "-c", f"{FIXTURES_DIR}/config/complete.json"])
        assert fetch.config == {
            "ena_api_username": "ENA_FAKE",
            "ena_api_password": "FAKE",
            "url_max_attempts": 10,
            "fire_endpoint": "fake_endpoint",
            "fire_ena_bucket": "fake_bucket",
            "fire_access_key_id": "",
            "fire_secret_access_key": "",
        }

    def test_config_override_partial_with_json(self):
        fetch = fetch_reads.FetchReads(argv=["-p", "ERP001736", "-c", f"{FIXTURES_DIR}/config/partial.json"])
        assert fetch.config == {
            "ena_api_username": "",
            "ena_api_password": "",
            "url_max_attempts": 8,
            "fire_endpoint": "fake_endpoint",
            "fire_ena_bucket": "fake_bucket",
            "fire_access_key_id": "",
            "fire_secret_access_key": "",
        }
