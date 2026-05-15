# nate-skills — Agent Instructions

This repository is a personal skills package for Claude Code, Codex, Gemini-CLI, and other
agent harnesses. Skills are reusable context bundles that agents load when working on a specific
codebase or domain.

`CLAUDE.md` is a symlink to this file (the cross-vendor convention).

## Repository layout

```
nate-skills/
├── AGENTS.md             ← this file (CLAUDE.md symlinks here)
├── README.md             ← human-facing index and usage guide
├── <skill-name>/
│   ├── SKILL.md          ← skill entry point (frontmatter + guide)
│   ├── <skill>.skill     ← compiled skill archive (ZIP)
│   └── references/
│       └── *.md          ← reference docs loaded on demand
└── local-agents/
    └── <skill-name>/     ← trimmed variants for smaller models
        ├── SKILL.md
        └── references/
```

## When a skill is added or updated

1. Put the skill in its own directory named after it (e.g., `my-skill/`).
2. Write `SKILL.md` with YAML frontmatter (`name`, `description`) followed by the skill body.
3. Put all reference documents in `my-skill/references/` — SKILL.md should reference them as
   `references/<file>.md`.
4. If trimmed variants for smaller models are warranted, mirror the structure under
   `local-agents/my-skill/` using `-small` (30–70B) or `-mini` (7–15B) suffixes in the
   filename frontmatter name.
5. Update the skills table in `README.md`.
6. Commit and push when done.

## Skill authoring conventions

- **SKILL.md frontmatter** must include `name` and `description`. The description is used by
  agent harnesses to match the skill to a task — make it specific and keyword-rich.
- **Reference files** should be self-contained and task-scoped. SKILL.md should tell the agent
  *which* reference to read *when*, not dump everything at once.
- **Compiled `.skill` archives** are ZIP bundles of the skill directory. Rebuild with the
  harness-specific build command if you restructure the directory.
- **Do not** add AL/Dynamics project boilerplate to `.gitignore` — keep it scoped to artifacts
  this repo actually produces (.DS_Store, build outputs, etc.).

## Error patterns to track

If you encounter a repeatable failure while using or building a skill (wrong path reference,
broken frontmatter, harness import error), document it here so future agents avoid it:

- *(none logged yet)*
