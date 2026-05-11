use std::path::Path;

use crate::{LibError, LibResult};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CommandSpec {
    program: String,
    args: Vec<String>,
}

impl CommandSpec {
    pub fn new(program: impl Into<String>, args: Vec<String>) -> LibResult<Self> {
        let program = program.into();
        validate_program(&program)?;
        for arg in &args {
            validate_arg(arg)?;
        }
        Ok(Self { program, args })
    }

    pub fn program(&self) -> &str {
        &self.program
    }

    pub fn args(&self) -> &[String] {
        &self.args
    }

    pub fn argv(&self) -> Vec<String> {
        let mut argv = Vec::with_capacity(self.args.len() + 1);
        argv.push(self.program.clone());
        argv.extend(self.args.clone());
        argv
    }
}

pub fn path_arg(path: &Path) -> LibResult<String> {
    let Some(value) = path.to_str() else {
        return Err(LibError::InvalidArguments(format!(
            "path is not valid UTF-8: {}",
            path.display()
        )));
    };
    validate_arg(value)?;
    Ok(value.to_owned())
}

fn validate_program(program: &str) -> LibResult<()> {
    if program.trim().is_empty() {
        return Err(LibError::InvalidArguments(
            "external tool program cannot be empty".to_owned(),
        ));
    }
    if has_shell_metachar(program) || program.contains('/') {
        return Err(LibError::InvalidArguments(format!(
            "external tool program must be a simple executable name: {program:?}"
        )));
    }
    Ok(())
}

fn validate_arg(arg: &str) -> LibResult<()> {
    if arg.contains('\0') {
        return Err(LibError::InvalidArguments(
            "external tool arguments cannot contain NUL bytes".to_owned(),
        ));
    }
    Ok(())
}

fn has_shell_metachar(value: &str) -> bool {
    value
        .chars()
        .any(|ch| matches!(ch, '|' | '&' | ';' | '<' | '>' | '`' | '$' | '\n' | '\r'))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_shell_programs() {
        assert!(CommandSpec::new("samtools;rm", vec![]).is_err());
        assert!(CommandSpec::new("/usr/bin/samtools", vec![]).is_err());
        assert!(CommandSpec::new("samtools", vec!["region;ok-as-arg".to_owned()]).is_ok());
    }
}
