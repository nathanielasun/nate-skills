# local-agents

Trimmed skill variants for local coding agents with smaller context windows.

## Model tiers

| Suffix | Parameter range | Notes |
|--------|----------------|-------|
| `-small` | 30–70B | Condensed references, same structure as frontier skill |
| `-mini` | 7–15B | Minimal core only — the most essential rules and patterns |

## Structure

Each skill mirrors its counterpart in the parent directory but with reduced content:

```
local-agents/
├── WarpX-cpp-small/
│   ├── SKILL.md
│   └── references/
├── WarpX-cpp-mini/
│   ├── SKILL.md
│   └── references/
├── WarpX-python-small/
│   ├── SKILL.md
│   └── references/
└── WarpX-python-mini/
    ├── SKILL.md
    └── references/
```

## Current variants

| Skill | Tier | Description |
|-------|------|-------------|
| [WarpX-cpp-small](WarpX-cpp-small/SKILL.md) | `-small` (30–70B) | WarpX C++ — condensed architecture, coding rules, collisions, build/test/PR |
| [WarpX-cpp-mini](WarpX-cpp-mini/SKILL.md) | `-mini` (7–15B) | WarpX C++ — minimal patterns, collisions, build/run |
| [WarpX-python-small](WarpX-python-small/SKILL.md) | `-small` (30–70B) | WarpX Python/PICMI — condensed input, collisions, output, runtime |
| [WarpX-python-mini](WarpX-python-mini/SKILL.md) | `-mini` (7–15B) | WarpX Python/PICMI — minimal PICMI basics, collisions, install/output |

## When to add a local-agent variant

Add a `-small` or `-mini` variant when the frontier skill's total token count (SKILL.md +
all references) exceeds the model's usable context for the task, or when running a local model
via Ollama, LM Studio, or similar.

Keep the frontmatter `name` field consistent — use the same base name with the suffix appended
(e.g., `warpx-cpp-small`, `warpx-cpp-mini`) so harnesses can discover variants automatically.
