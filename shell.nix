{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = with pkgs; [
    # Rust toolchain
    cargo
    rustc
    rustfmt
    clippy
    rust-analyzer
    
    # Build dependencies (common ones, adjust as needed)
    pkg-config
    
    # Additional tools
    cargo-watch
    cargo-edit
  ];

  # Environment variables
  RUST_BACKTRACE = 1;
}
