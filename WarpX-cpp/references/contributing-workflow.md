# Contributing workflow: PRs, documentation, review expectations

This reference covers the WarpX contribution workflow end-to-end — git setup, branching, PR preparation, documentation requirements, the interaction with WarpX's own in-repo skills, and (importantly) the community's explicit expectations around LLM-generated code. Read this before opening a PR.

## Table of contents

1. [The 30-second summary](#the-30-second-summary)
2. [Git setup and forking](#git-setup-and-forking)
3. [Branching and the development workflow](#branching-and-the-development-workflow)
4. [Writing the commit messages](#writing-the-commit-messages)
5. [Opening the pull request](#opening-the-pull-request)
6. [PR sizing: small, focused, separable](#pr-sizing-small-focused-separable)
7. [Documentation requirements](#documentation-requirements)
8. [Building and previewing documentation locally](#building-and-previewing-documentation-locally)
9. [Interaction with WarpX's in-repo skills](#interaction-with-warpxs-in-repo-skills)
10. [Updating AGENTS.md / CLAUDE.md](#updating-agentsmd--claudemd)
11. [LLM-generated code: explicit community expectations](#llm-generated-code-explicit-community-expectations)
12. [Code review: what to expect](#code-review-what-to-expect)
13. [After the PR is merged](#after-the-pr-is-merged)

## The 30-second summary

For any change beyond a typo fix:

1. Fork `BLAST-WarpX/warpx` on GitHub.
2. Create a feature branch off `development` in your fork.
3. Make changes in the smallest reasonable scope.
4. Add a test (in `Tests/Tests/`) and documentation (in `Docs/source/`) in the *same* PR as the feature.
5. Run `clang-tidy` and a local build (at minimum 2D, ideally 3D + GPU).
6. Open a PR against `BLAST-WarpX/warpx:development`.
7. Wait for CI; respond to review feedback; iterate.
8. If LLM-assisted: read the code carefully yourself before requesting review.

The full elaboration follows.

## Git setup and forking

One-time setup for your GitHub account:

1. **Real name and affiliation**: set on your GitHub profile at `https://github.com/settings/profile`. Reviewers want to know who they're talking to.
2. **Verified email**: add your professional email at `https://github.com/settings/emails`.
3. **Local git config** matching the GitHub account:
   ```bash
   git config --global user.name "FIRSTNAME LASTNAME"
   git config --global user.email EMAIL@EXAMPLE.com
   ```
4. **SSH key** at `https://github.com/settings/keys`.

Fork the repository at `https://github.com/BLAST-WarpX/warpx` (button on the GitHub UI). This creates a personal copy you have write access to.

Clone setup (one-time):

```bash
# Clone the mainline (read-only for you)
git clone git@github.com:BLAST-WarpX/warpx.git
cd warpx

# Rename "origin" to "mainline" to make it clear
git remote rename origin mainline

# Add your fork as a separate remote
git remote add myGithubUsername git@github.com:myGithubUsername/warpx.git

# Verify
git remote -v
# should show:
#   mainline  git@github.com:BLAST-WarpX/warpx.git (fetch/push)
#   myGithubUsername  git@github.com:myGithubUsername/warpx.git (fetch/push)
```

This setup distinguishes "pull from the project" (`git pull mainline development`) from "push my work" (`git push myGithubUsername feature-branch`).

## Branching and the development workflow

All new work starts from an up-to-date `development` branch:

```bash
git checkout development
git pull mainline development

# Create a feature branch with a descriptive name
git checkout -b add-new-mcc-process
```

Branch naming: descriptive, hyphenated, related to the change. Examples:
- `fix-mcc-energy-conservation`
- `add-coulomb-collisions-rz`
- `refactor-multifab-register-access`
- `docs-clarify-pusher-options`

Avoid branch names like `dev`, `mywork`, `feature1` — they don't help reviewers or your future self.

When the `development` branch moves while you're working (common — the project is active), sync your branch:

```bash
git checkout add-new-mcc-process
git fetch mainline
git merge mainline/development
# or, if you prefer rebasing:
# git rebase mainline/development
```

Resolve any conflicts. WarpX maintainers generally prefer merges over rebases for shared branches, but for your own feature branch either is fine — pick what you're comfortable with.

## Writing the commit messages

The first line of a commit message is the title: short, imperative, ~40-60 characters. Then a blank line, then a body if needed (paragraphs explaining *why*, not *what* — the diff shows what).

```
Add charge-exchange MCC process for noble-gas ions

The existing MCC infrastructure handled electron-neutral collisions
but did not have a clean path for ion-neutral processes. This adds
ChargeExchange.H following the pattern of ImpactIonization.H, wires
it into BackgroundMCCCollision::doCollisions, and adds a 2D test
case using krypton (the canonical low-temperature-plasma reference
gas).

The cross-section file follows the LXCat convention; the cross
section for Kr+ + Kr was taken from Phelps (1991).
```

Don't:
- Write commit messages like `wip`, `fix typo`, `changes`. These obscure history.
- Include a long bullet list — the diff shows the bullets.
- Reference issue numbers in the title line only — put them in the PR description instead.

It's fine to have many small commits during development. Before opening a PR, consider `git rebase -i` to squash trivial commits and produce a clean history. But this is optional — reviewers often look at the PR as a whole, not commit-by-commit.

## Opening the pull request

After pushing to your fork:

```bash
git push -u myGithubUsername add-new-mcc-process
```

GitHub shows a banner inviting you to open a PR. Click "Compare & pull request."

In the PR description, include:

- **What the change does** (one paragraph, suitable for the release notes).
- **Why it's needed** (the motivation, the issue/discussion it resolves).
- **How you tested it** (which tests pass, any new tests added).
- **Anything reviewers should pay attention to** (tricky parts, areas you're uncertain about, follow-up work).
- **Screenshots or plots** if the change is user-visible (new diagnostic outputs, new physics validation).
- **Reference links** to related issues, discussions, papers, prior PRs.

PR title format: a short descriptive line. No specific prefix required, but `[WIP]` (Work-In-Progress) is the conventional way to mark PRs that are not yet ready for review:

```
[WIP] Add charge-exchange MCC process for noble-gas ions
```

Remove the `[WIP]` once the PR is ready for review. Reviewers will (rightly) ignore WIP PRs unless explicitly tagged for early feedback.

The PR target branch is always **`BLAST-WarpX/warpx:development`**, not `main` or `release-*`. Releases are cut from `development` periodically.

## PR sizing: small, focused, separable

The official guidance: "Please DO NOT write large pull requests, as they are very difficult and time-consuming to review. As much as possible, split them into small, targeted PRs."

Examples of good splits:

- Adding a new physics module is *three* PRs: (1) infrastructure refactor that the module depends on, (2) the module itself with one minimal test, (3) extended tests and edge-case handling.
- Adding a new MCC process plus a documentation cleanup of the existing MCC docs is *two* PRs: process first, docs separately.
- Fixing a bug and refactoring the surrounding code is *two* PRs: minimal bug fix first, refactor afterward.
- Style changes are *always* a separate PR from substantive changes. "It is recommended that style changes are not included in the PR where new code is added. This is to avoid any errors that may be introduced in a PR just to do style change."

If you're not sure how to split, open an issue first describing the planned work and ping the maintainers — they'll often suggest a decomposition. This is much cheaper than discovering the wrong split during review.

PRs with more than ~500-1000 lines of code changed are hard to review carefully. If you have more, split into a stack of dependent PRs or merge a foundation PR first and then build on it.

## Documentation requirements

Every new feature must include documentation in the same PR. The docs are in `Docs/source/` and use Sphinx (reStructuredText, `.rst`).

The documentation locations to know:

| Change type | Documentation location |
|---|---|
| New input file parameter | `Docs/source/usage/parameters.rst` |
| New physics module / theory | `Docs/source/theory/<topic>.rst` (e.g., `multiphysics/collisions.rst`) |
| New algorithm in an existing module | `Docs/source/theory/models_algorithms/<topic>.rst` |
| New numerical method | Possibly its own section under `theory/` |
| New developer-facing implementation detail | `Docs/source/developers/<topic>.rst` |
| New how-to workflow | `Docs/source/usage/workflows/<name>.rst` |
| New example | `Docs/source/usage/examples/<name>/` with a README and input file |
| New highlight (paper using WarpX) | `Docs/source/highlights.rst` |
| New tutorial | `Docs/source/tutorials/<name>.rst` |

For runtime parameters, document them in `parameters.rst` with the syntax:

```rst
* ``my_module.my_parameter`` (`bool`) optional (default: ``0``)
    Brief description of what the parameter does, with edge cases and dependencies.

* ``my_module.other_parameter`` (`float`) optional (default: ``1.0``)
    Description here.
```

For C++ APIs that are user-facing or worth highlighting, the Doxygen blocks in the source code feed the auto-generated reference docs. Use this format:

```cpp
/**
 * @brief One-line summary of what the function does.
 *
 * Longer description if needed, including any non-obvious behavior or
 * preconditions.
 *
 * @param[in] x  what x is, units, valid range
 * @param[out] result  what gets returned
 * @return what the function returns and what it means
 */
amrex::Real my_function (amrex::Real x);
```

## Building and previewing documentation locally

To verify your documentation renders correctly:

```bash
cd Docs
make html
```

This requires Sphinx, the WarpX-specific Sphinx extensions, and the RTD theme. Install with:

```bash
pip install -r Docs/requirements.txt
```

The output appears in `Docs/build/html/index.html`. Open in a browser to navigate the rendered docs.

Common issues:
- **Cross-references that don't resolve**: `:ref:` and `:doc:` directives need correct target labels. Check the warnings in `make html` output.
- **Math rendering**: math is rendered via MathJax in the browser. Test in a real browser, not just by checking the HTML source.
- **Indentation in RST**: RST is sensitive to indentation in directives. Look at adjacent existing content for the indentation level.

For one-off doc edits, you can also use the GitHub web UI's preview, but `make html` is more reliable.

## Interaction with WarpX's in-repo skills

WarpX maintains its own LLM-facing tooling at `.claude/skills/`. The two skills as of the early 2026 development branch:

### `/warpx-answer-user-question`

A skill for drafting responses to user questions posted in WarpX GitHub issues, discussions, or emails. The workflow:

```
/warpx-answer-user-question https://github.com/BLAST-WarpX/warpx/discussions/1234
```

The skill fetches the question, categorizes it (installation, input parameters, physics, etc.), searches the source code and docs, drafts a response in the style of an experienced developer, and presents the draft for human review.

If your task is "help me reply to a user question," prefer this in-repo skill to writing the response from scratch.

### `/warpx-new-paper-highlight`

Adds a new paper to `Docs/source/highlights.rst`. Usage:

```
/warpx-new-paper-highlight https://doi.org/10.1103/PhysRevLett.133.045002
```

Fetches paper metadata, picks the appropriate highlights section (Plasma-Based Acceleration, HPC and Numerics, etc.), formats the entry in the RST style of the file, and optionally creates a branch and PR.

If your task is "add a paper to highlights," use this in-repo skill rather than editing `highlights.rst` manually.

### When to use the in-repo skills vs this skill

- Use the **in-repo skills** for the workflows they're designed for (user-question replies, paper highlights). They have version-matched context and produce output in the project's exact style.
- Use **this skill** (`warpx-cpp`) for: writing new C++ code, debugging, navigating the codebase, refactoring, adding tests for C++ features, understanding the architecture.
- The skills are complementary, not exclusive. A larger task may use both — e.g., add a paper highlight (use `/warpx-new-paper-highlight`) and add the C++ feature the paper describes (use this skill).

## Updating AGENTS.md / CLAUDE.md

`AGENTS.md` at the repo root contains the LLM-facing project instructions. It's kept under 300 lines for context-window efficiency. `CLAUDE.md` is a symlink to `AGENTS.md` (both names are read by different LLM clients).

When to update `AGENTS.md`:
- A new build command or option becomes essential.
- A new conventional pattern emerges that LLMs should follow.
- An existing instruction becomes outdated or wrong.
- A new in-repo skill or workflow is added.

How to update:
- Keep the file concise. Aim for instructions that change LLM behavior in observable, beneficial ways — not just nice-to-haves.
- Use imperative voice. "Always do X. Never do Y. Prefer Z over W."
- Don't duplicate what's in the documentation. Point to it.
- Keep under 300 lines. If you're approaching that limit, remove something less important.

`AGENTS.md` is part of the repository, so it goes through the same PR process as code changes. PRs touching `AGENTS.md` should explain why the change improves LLM-assisted development.

## LLM-generated code: explicit community expectations

The WarpX community has explicit, documented expectations for code generated with LLM assistance. From the official `how_to_develop_with_llms.html` page:

> The WarpX community thus urges you to perform a careful, manual review of all LLM-generated code and documentation before asking for a review of your pull-request. This is important, otherwise you risk to waste the valuable time of our most proficient developers that will need to review your LLM-generated code.
>
> (Be considerate that WarpX developers can prompt an LLM just as efficiently as you can. Your critical thinking skill to make sense of the LLM-generated code and make it sensible for review and maintainable for the long term is what is needed!)

What this means in practice:

1. **You read every line of LLM-generated code before submitting it.** Not skim — read, understand, and be prepared to defend each design decision in review.
2. **You verify against the actual source, not just the LLM's plausible-sounding claim.** LLMs hallucinate APIs, miss physics constraints, and produce code that compiles but is subtly wrong. The skill files in `warpx-cpp` (this skill) include explicit caveats about API stability for exactly this reason.
3. **You write the test before requesting review.** Tests are a sanity check on both you and the LLM. A test that passes for a wrong reason is a worse outcome than no test at all — verify the test exercises what you think it does.
4. **You attribute the LLM assistance.** It's not required to flag every line as LLM-generated, but if a major part of a PR was produced with substantial LLM help, mention it in the PR description so reviewers can calibrate.
5. **You don't open low-quality PRs.** A PR that the maintainer can tell was minimally reviewed by the contributor wastes everyone's time. The bar is "would you have opened this PR if you'd written it yourself by hand?" If not, don't open it.

This skill (warpx-cpp) is built around supporting these expectations — pointing to source for verification, recommending precedent-based development, flagging API stability concerns. Use it that way.

## Code review: what to expect

After opening a PR, expect:

1. **CI runs automatically.** Multiple jobs across geometries, compute backends, and operating systems. Tests must pass. clang-tidy must not introduce new warnings. Documentation must build.

2. **A maintainer assigns themselves or pings the right person.** Response time varies — for a focused PR, days to a week is typical. Large or unusual PRs may take longer.

3. **Review comments arrive in batches.** Most reviewers do a single pass with multiple suggestions rather than commenting line-by-line in real time. Read all comments before responding.

4. **You address each comment.** Either:
   - Make the requested change, push a new commit, and mark the conversation resolved.
   - Reply explaining why you think the existing code is correct, and discuss with the reviewer.
   
   Don't silently ignore comments. If you disagree, say so explicitly.

5. **Multiple review rounds are normal.** Even small PRs often go through 2-3 cycles. Don't take this personally — it's how the project maintains quality.

6. **Approvals from at least one maintainer** (and sometimes from a domain expert if the PR touches specialized physics) are required before merging.

7. **Squash-and-merge or merge commit** depending on the PR. The maintainer usually chooses. Don't worry about cleaning history beyond what's reasonable; the merge will be tidy regardless.

Things that speed up review:
- A clear PR description.
- Small PRs.
- Tests that exercise the change.
- Documentation included.
- Responses to comments within a week.
- Visible care in the code (no dead code, no commented-out experiments, no debug `Print` statements left in).

Things that slow it down (or cause rejection):
- LLM-generated code that the contributor obviously hasn't read carefully.
- Style violations (missing `m_` prefix, raw `0.0` literals, `using namespace amrex`).
- A test that doesn't actually test the new code.
- Mixing style changes with substantive changes.
- "Just merge this, we'll fix it later" energy. Reviewers won't.
- PRs that touch large areas of code without prior discussion.

## After the PR is merged

When your PR merges:

1. **Delete the branch on your fork** (GitHub offers a button after merge). Keeps your repository tidy.
2. **Sync your local `development`** with the mainline:
   ```bash
   git checkout development
   git pull mainline development
   git branch -d add-new-mcc-process   # delete local branch
   ```
3. **Watch for follow-up issues.** If your feature has rough edges or limitations you flagged in the PR description, file follow-up issues so they don't get lost.
4. **Update related documentation** if the merge revealed gaps. Often a doc PR follows a feature PR by a few days.
5. **If your contribution is associated with a paper**, consider opening a `/warpx-new-paper-highlight` PR once the paper is published.

That's the workflow. Each step is documented to lower the activation energy for new contributors and to keep the bar high enough that maintainers' time isn't wasted. Respect both ends of that bargain.

## A note on emergency or large-scale changes

For changes that don't fit the normal workflow:

- **Critical bug fixes affecting current users**: open an issue first, ping maintainers, then a focused fix PR. These can move faster than feature PRs.
- **Large refactors touching many files**: open a discussion at `https://github.com/BLAST-WarpX/warpx/discussions` first to socialize the design. Don't surprise maintainers with a 5,000-line refactor in a cold PR.
- **API breakage**: requires explicit maintainer sign-off in a discussion or issue. Communicate in advance.
- **Dependency version bumps** (AMReX, PICSAR, openPMD, pyAMReX): these are typically done by maintainers on a regular cadence. If you need a feature from a newer dependency, ask whether the bump can be moved up.

When in doubt, open a discussion or issue before opening a PR. The cost of a clarifying conversation is much lower than the cost of revising a large PR that went the wrong direction.
