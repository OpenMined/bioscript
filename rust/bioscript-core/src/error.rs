use std::{error::Error, fmt};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RuntimeError {
    Monty(String),
    Unsupported(String),
    InvalidArguments(String),
    Io(String),
}

impl fmt::Display for RuntimeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Monty(msg)
            | Self::Unsupported(msg)
            | Self::InvalidArguments(msg)
            | Self::Io(msg) => f.write_str(msg),
        }
    }
}

impl Error for RuntimeError {}
