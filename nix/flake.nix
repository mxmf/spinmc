{
  description = "Rust devshell for Iced";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
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
          pname = "MC-Curie";
          version = "0.1.0";

          src = ./..;

          cargoLock = {
            lockFile = ../Cargo.lock;
          };

          buildInputs = [ ];

          nativeBuildInputs = [ pkgs.pkg-config ];

          meta = with pkgs.lib; {
            description = "mc_curie";
            license = licenses.mit;
            maintainers = with maintainers; [ mxmf ];
          };
        };

        devShells.default = mkShell rec {
          buildInputs = [
            rust-analyzer
            rust-bin.stable.latest.default
            cargo-flamegraph
            cargo-expand
            cargo-show-asm

            hdf5
            pkg-config


            python3
            python3Packages.matplotlib
            python3Packages.scipy
            python3Packages.pyqt5
            qt5.qtwayland
            libsForQt5.qt5.qtbase


            # pkg-config
          ] ++ lib.optionals (stdenv.isLinux) [
          ];

          MPLBACKEND = "QtAgg";


          LD_LIBRARY_PATH = lib.makeLibraryPath buildInputs;


        };
      }
    );
}
