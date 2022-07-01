#
# Copyright 2020 Adobe. All rights reserved.
# This file is licensed to you under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. You may obtain a copy
# of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under
# the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR REPRESENTATIONS
# OF ANY KIND, either express or implied. See the License for the specific language
# governing permissions and limitations under the License.
#
if(TARGET nlohmann_json::nlohmann_json)
    return()
endif()

message(STATUS "Third-party: creating target 'nlohmann_json::nlohmann_json'")

# nlohmann_json is a big repo for a single header, so we just download the release archive
include(FetchContent)
FetchContent_Declare(
    nlohmann_json
    URL "https://github.com/nlohmann/json/releases/download/v3.10.5/json.tar.xz"
    URL_HASH SHA256=344be97b757a36c5b180f1c8162f6c5f6ebd760b117f6e64b77866e97b217280
)
FetchContent_MakeAvailable(nlohmann_json)