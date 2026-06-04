#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObjectKind {
    Module,
    AlignmentFile,
    AlignedSegment,
    Fasta,
    FastaRecord,
    VariantFile,
    VariantRecord,
}

#[derive(Debug, Clone, PartialEq)]
pub enum LibValue {
    None,
    Bool(bool),
    Int(i64),
    String(String),
    List(Vec<LibValue>),
    Object(ObjectKind),
}
