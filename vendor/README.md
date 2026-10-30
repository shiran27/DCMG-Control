# MATLAB SDP solver

The local `sedumi/` folder is the unmodified official SeDuMi 1.3.8 release:
https://github.com/sqlp/sedumi/releases/tag/v1.3.8

It includes Apple Silicon MATLAB binaries and its GPL license. It is a local
dependency, ignored by Git. The experiment prefers a working MOSEK installation
and falls back to SeDuMi or SDPT3. It never modifies the saved MATLAB path.

To reproduce on another computer, download and extract `sedumi.tgz` from that
release into this folder, or install a supported SDP solver on the MATLAB path.
YALMIP is also required. Use `install_sedumi -rebuild` inside the extracted
folder only if the bundled binaries do not support your platform.

Downloaded archive SHA-256:
`9502d62812f8b99324b06e19a95ccb45753efff65c251608b11538dad9d08300`.
