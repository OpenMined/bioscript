# expo-bioscript

Expo module wrapper for the BioScript runtime.

## Current

BioScript currently lives in this repo as a Rust-first runtime for secure genomic analysis with Pythonic syntax via Monty.

Right now:
- `rust/bioscript` contains the runtime logic
- `monty/` is the interpreter/runtime dependency
- `expo-bioscript/` contains a working wrapper scaffold with Android and iOS packaging scripts
- BioScript is still evolving, so the mobile-facing API is not treated as stable

The current goal is to expose a narrow Expo-native interface to BioScript without exposing generic Monty execution directly.

## Near-Term

The expected implementation is:
- keep Monty internal
- add a thin Rust FFI layer for mobile-safe entrypoints
- expose a small Expo API such as:
  - `runFile(...)`
  - `runCode(...)`
  - `prepareIndexes(...)`
- enforce resource limits, restricted I/O, and blocked OS/network access through the BioScript runtime

This `expo-bioscript` folder is the home for the Expo wrapper layer:
- `ios/`
- `android/`
- `src/`
- `scripts/`

Current wrapper status:
- `runFile(...)` is implemented end-to-end through the Rust FFI layer
- Android native packaging builds successfully
- iOS native packaging builds successfully
- Apple mobile targets currently disable HTS-backed CRAM/BAM indexing and lookup paths

## First API

The first Expo-facing API should be `runFile(...)`, not `runCode(...)`.

Reason:
- it fits the current BioScript runtime better
- it keeps execution rooted in explicit files instead of arbitrary code strings
- it is the safer first bridge for mobile

Proposed request shape:

```ts
type RunFileRequest = {
  scriptPath: string;
  root?: string;
  inputFile?: string;
  outputFile?: string;
  participantId?: string;
  traceReportPath?: string;
  timingReportPath?: string;
  inputFormat?: 'auto' | 'text' | 'zip' | 'vcf' | 'cram';
  inputIndex?: string;
  referenceFile?: string;
  referenceIndex?: string;
  autoIndex?: boolean;
  cacheDir?: string;
  maxDurationMs?: number;
  maxMemoryBytes?: number;
  maxAllocations?: number;
  maxRecursionDepth?: number;
};
```

Proposed response shape:

```ts
type RunFileResult = {
  ok: true;
};
```

For the first implementation, errors should be surfaced as native/module exceptions rather than encoded into the success payload.

## Long-Term

If the API stabilizes, this will likely evolve into a cleaner split:

- `bioscript-core`
  - runtime, domain logic, Monty integration
- `bioscript-ffi`
  - stable native/mobile boundary
- `expo-bioscript`
  - Expo module wrapper

At that point, `expo-bioscript` may move into its own repo and depend on BioScript through a more stable versioned interface.
