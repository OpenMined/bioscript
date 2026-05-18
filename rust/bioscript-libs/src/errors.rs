use std::fmt;

pub type LibResult<T> = Result<T, LibError>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LibError {
    UnknownModule(String),
    UnsupportedMode {
        module: &'static str,
        object: &'static str,
        mode: String,
    },
    UnsupportedFeature {
        module: &'static str,
        feature: &'static str,
    },
    InvalidArguments(String),
}

impl LibError {
    pub fn unsupported_feature(module: &'static str, feature: &'static str) -> Self {
        Self::UnsupportedFeature { module, feature }
    }
}

impl fmt::Display for LibError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownModule(name) => write!(f, "unknown bioscript library module: {name}"),
            Self::UnsupportedMode {
                module,
                object,
                mode,
            } => write!(f, "{module}.{object} does not support mode {mode:?}"),
            Self::UnsupportedFeature { module, feature } => {
                write!(f, "{module} does not support {feature}")
            }
            Self::InvalidArguments(message) => write!(f, "{message}"),
        }
    }
}

impl std::error::Error for LibError {}
