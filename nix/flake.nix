{
  description = "Rust devshell for Iced";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-26.05";
    rust-overlay.url = "github:oxalica/rust-overlay";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { nixpkgs, rust-overlay, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
      in
      with pkgs;
      {
        packages.default = pkgs.rustPlatform.buildRustPackage {
          pname = "spinmc";
          version = "0.1.2";

          src = ./..;

          cargoLock = {
            lockFile = ../Cargo.lock;
          };

          buildInputs = [ ];

          nativeBuildInputs = [ pkgs.pkg-config ];

          meta = with pkgs.lib; {
            description = "spinmc";
            license = licenses.mit;
            maintainers = with maintainers; [ mxmf ];
          };
        };

        devShells.default = mkShell (
          rec {
            buildInputs = [
              rust-analyzer
              rust-bin.stable.latest.default
              cargo-flamegraph
              cargo-expand
              cargo-show-asm

              hdf5
              cmake
              pkg-config

              uv
              maturin
            ] ++ lib.optionals (stdenv.isLinux) [
              qt5.qtwayland
              libsForQt5.qtbase
              python3
              python3Packages.matplotlib
              python3Packages.pyqt5
            ];

            CMAKE_POLICY_VERSION_MINIMUM = "3.5";

            MPLBACKEND = lib.optionalString stdenv.isLinux "QtAgg";

            LD_LIBRARY_PATH = lib.makeLibraryPath buildInputs;

            shellHook = ''
              if [ ! -d ".venv" ]; then
                uv venv .venv
              fi
              source .venv/bin/activate
              echo "uv pip env ready"
            '';
          }
          // lib.optionalAttrs stdenv.isDarwin {
            # nix's cc-wrapper targets aarch64-apple-darwin, but chemfiles-sys's
            # cmake build passes --target arm64-apple-macosx; align them here.
            CMAKE_OSX_ARCHITECTURES = "arm64";
            MACOSX_DEPLOYMENT_TARGET = "14.0";
          }
        );
      }
    );
}
