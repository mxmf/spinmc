fn main() {
    if std::env::var_os("CARGO_CFG_WINDOWS").is_some() {
        // chemfiles-sys references GetUserNameA but does not currently link the
        // Windows system library that provides it.
        println!("cargo:rustc-link-lib=Advapi32");
    }
}
