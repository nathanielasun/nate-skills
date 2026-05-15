# nate-skills

Personal skills package for Claude Code, Codex, Gemini-CLI, and other agent harnesses.
Skills are reusable context bundles that agents load when working on a specific codebase or domain.

## Available skills

| Skill | Description | Models |
|-------|-------------|--------|
| [warpx-cpp](WarpX-cpp/SKILL.md) | WarpX C++ codebase — particle-in-cell physics, AMReX architecture, MCC/DSMC collisions, GPU portability, PR workflow | Frontier |
| [warpx-python](WarpX-python/SKILL.md) | WarpX Python/PICMI interface — input scripts, species/injection, diagnostics, runtime callbacks, post-processing | Frontier |

Local-model variants (smaller context, trimmed references) live in [local-agents/](local-agents/) — see that directory for the full list of `-small` (30–70B) and `-mini` (7–15B) variants.

---

## Importing a skill

### Claude Code
```bash
# From a local checkout of this repo:
claude skill add /path/to/nate-skills/WarpX-cpp/warpx-cpp.skill

# Or reference the SKILL.md directly in a project's .claude/skills/ directory
cp -r /path/to/nate-skills/WarpX-cpp .claude/skills/warpx-cpp
```

### Codex / OpenAI Agents
Point the agent's context config at the `SKILL.md` file — Codex reads the YAML frontmatter
`description` field to decide when to activate the skill.

### Gemini-CLI
Use `gemini context add <path/to/SKILL.md>` or include the skill directory in your project's
`.gemini/context/` folder.

### Generic / other harnesses
Load `SKILL.md` as a system-prompt prefix. The frontmatter `description` tells the agent when
to activate; the body tells it what to do.

---

## Adding a new skill

1. **Create a directory** named after the skill: `my-skill/`
2. **Write `SKILL.md`** — YAML frontmatter at the top, then the skill body:
   ```yaml
   ---
   name: my-skill
   description: Trigger description — what tasks/keywords should activate this skill.
   ---
   ```
3. **Put reference docs in `my-skill/references/`** and reference them from SKILL.md as
   `references/<file>.md`. Keep each reference file task-scoped so the agent only reads
   what it needs.
4. **Add local-model variants** (optional) under `local-agents/my-skill/` with `-small`
   (30–70B) or `-mini` (7–15B) suffixes in the frontmatter `name`.
5. **Update the table above** with a one-line description.
6. **Compile the `.skill` archive** if your harness requires it (ZIP of the skill directory).
7. Commit and push.

See [AGENTS.md](AGENTS.md) for the full authoring conventions and error-pattern log.
