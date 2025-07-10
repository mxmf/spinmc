{
  description = "Iced example";

  inputs = {
    flake-parts.url = "github:hercules-ci/flake-parts";
    cargo2nix.url = "github:cargo2nix/cargo2nix/release-0.12";
    rust-overlay.url = "github:oxalica/rust-overlay";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = inputs@{ self, nixpkgs, flake-parts, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = nixpkgs.lib.systems.flakeExposed;
      perSystem = { self', pkgs, system, ... }:
        let
          rustVersion = "1.88.0";
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ inputs.cargo2nix.overlays.default (import inputs.rust-overlay) ];
          };
          runtimeDeps = with pkgs; [
            wayland
            libxkbcommon
            libGL
          ];

          rustPkgs = pkgs.rustBuilder.makePackageSet {
            inherit rustVersion;
            packageFun = import ./Cargo.nix;
            packageOverrides = pkgs: pkgs.rustBuilder.overrides.all ++ [
              (pkgs.rustBuilder.rustLib.makeOverride {
                name = "counter";
                overrideAttrs = drv: {
                  propagatedNativeBuildInputs = drv.propagatedNativeBuildInputs or [ ] ++ (with pkgs; [
                    pkg-config
                    makeWrapper
                  ] ++ runtimeDeps);
                  postFixup = ''
                    patchelf --set-rpath ${pkgs.lib.makeLibraryPath runtimeDeps} $bin/bin/counter
                  '';
                };
              })
            ];
          };
        in
        {
          packages = rec {
            counter = (rustPkgs.workspace.counter { }).bin;
            default = counter;
          };
          devShells.default = pkgs.mkShell rec {
            buildInputs = with pkgs; [
              pkg-config
            ] ++ runtimeDeps ++ [
              rust-analyzer-unwrapped
              (rust-bin.stable.${rustVersion}.default.override { extensions = [ "rust-src" ]; })

              glib
              atk
              gtk3
              # gtkd
              clang
              protobuf
              xdo
              xdot
              xdotool
              libayatana-appindicator

            ];
            LD_LIBRARY_PATH = "${nixpkgs.lib.makeLibraryPath buildInputs}";
            LIBCLANG_PATH = "${pkgs.libclang.lib}/lib";

          };
        };
    };
}
