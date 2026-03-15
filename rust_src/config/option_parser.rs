use std::collections::HashMap;

pub type OptionMap = HashMap<String, String>;

pub fn parse_options(args: &[String]) -> Result<OptionMap, String> {
    let mut options = OptionMap::new();
    let mut i = 1;
    while i < args.len() {
        let key = args[i].clone();
        if i + 1 < args.len() && !args[i + 1].starts_with('-') {
            options.insert(key, args[i + 1].clone());
            i += 2;
        } else {
            options.insert(key, "true".to_string());
            i += 1;
        }
    }
    Ok(options)
}
