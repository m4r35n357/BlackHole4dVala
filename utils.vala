/*
Copyright (c) 2014, 2015, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using Json;

namespace Simulations {

    /**
     * Read JSON data from stdin
     */
    private static string fromStdin () {
        var input = new StringBuilder();
        var buffer = new char[1024];
        while (!stdin.eof()) {
            var read_chunk = stdin.gets(buffer);
            if (read_chunk != null) {
                input.append(read_chunk);
            }
        }
        return input.str;
    }

    /**
     * Parse JSON initial conditions data
     */
    public static Json.Object getJson () {
        unowned Json.Object obj;
        Json.Parser parser = new Json.Parser();
        try {
            parser.load_from_data(fromStdin());
            obj = parser.get_root().get_object();
        } catch (Error e) {
            stderr.printf("Unable to parse the data file: %s\n", e.message);
            return new Json.Object();
        }
        return obj;
    }

    public static int main (string[] args) {
        var arg0 = args[0].split("/");
        switch (arg0[arg0.length - 1]) {  // basename
            case "bh3d":
                KerrGeodesic.fromJson().solve();
                break;
            case "newton":
                Newton.fromJson().solve();
                break;
            case "nbody3d":
                NBody.fromJson().solve();
                break;
            default:  // for debugging "utils" binary
                switch (args[1]) {  // command line argument
                    case "bh3d":
                        KerrGeodesic.fromJson().solve();
                        break;
                    case "newton":
                        Newton.fromJson().solve();
                        break;
                    case "nbody3d":
                        NBody.fromJson().solve();
                        break;
                    default:
                        stdout.printf("Please specify an executable by name or by program argument");
                        break;
                }
                break;
        }
        return 0;
    }
}

