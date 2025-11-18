# AI Agent Guide for Chemparser Project

## What Are AI Agents?

AI agents are specialized sub-processes that can work autonomously on specific tasks. Think of them as junior developers you can assign tasks to - they'll work independently and report back with results.

## Available Agent Types

### 1. **general-purpose**
Best for: Complex multi-step tasks, research, code changes
- Full access to all tools (Read, Write, Edit, Bash, etc.)
- Can search code, make changes, run tests
- Use for: Feature implementation, bug fixes, refactoring

### 2. **Explore**
Best for: Understanding codebases quickly
- Specialized in finding files and searching code
- Faster than general-purpose for exploration
- Use for: "Where is X defined?", "How does Y work?"

### 3. **Plan**
Best for: Planning implementation approaches
- Creates step-by-step plans
- Doesn't execute code
- Use for: Breaking down complex features

## How to Use Agents

### Basic Pattern

```
I want you to [task description].
Use an agent to work on this autonomously.
```

### Best Practices

1. **Be Specific**: Give agents clear, detailed instructions
2. **Provide Context**: Tell them which files are relevant
3. **Set Expectations**: Explain what "done" looks like
4. **Launch Multiple**: Run agents in parallel for independent tasks

## Current Project Agent Configuration

### Auto-Approval Status
✅ ALL TOOLS AUTO-APPROVED via `.claude/settings.json`

This means agents can:
- Read/write/edit any file
- Run any bash command
- Search code
- Make git commits (if you ask)
- Launch other agents

### Safety Note
Auto-approval is safe for this project because:
- Agents are given specific tasks
- You review their changes after completion
- Git tracks all changes
- Agents won't do anything you didn't ask for

## Effective Agent Usage for Chemparser

### Example 1: Feature Implementation

**Bad:**
```
Add image size controls
```

**Good:**
```
Add image size controls to the Chrome extension.

Requirements from Todolist.md lines 4-6:
- Up/down arrows in bottom left of each image
- Developer option to save per-image size
- Developer option to save per-molecule size (all instances)

Files to modify:
- chem-extension/content.js (add UI + logic)
- chem-extension/popup.html (add options)
- chem-extension/popup.js (handle settings)

Create test file to verify it works.
```

### Example 2: Parallel Agent Execution

You can launch multiple agents at once:

```
Launch 3 agents in parallel:
1. Agent 1: Implement PubChem integration
2. Agent 2: Create 10 UI variations for popup
3. Agent 3: Optimize cache deduplication

Each agent should work independently and report results.
```

### Example 3: Research Tasks

```
Use an Explore agent to:
1. Find all places where SVG colors are modified
2. Document the color transformation logic
3. Identify any inconsistencies

Report findings in a markdown file.
```

## Common Agent Tasks for This Project

### Extension Development
```
Agent Task: Add [feature] to Chrome extension
Context: Read .claude.md and chem-extension/content.js
Requirements: [specific details from Todolist.md]
Deliverables: Modified code + test file
```

### Server Development
```
Agent Task: Implement [server feature]
Context: Read CLAUDE_CONTEXT.md
Files: [specific server files]
Test: Create curl commands to verify
```

### UI Design
```
Agent Task: Create UI variation #[N] for popup
Style: [description of desired look]
Requirements: Keep all existing functionality
Files: chem-extension/popup.html, styles.css
```

### Bug Fixes
```
Agent Task: Fix [bug description]
Reproduction: [how to reproduce]
Expected: [what should happen]
Current: [what actually happens]
```

## Monitoring Agent Progress

Agents will show their work in real-time:
- Tool calls they're making
- Files they're reading/editing
- Commands they're running

You'll see output like:
```
📖 Reading chem-extension/content.js...
✏️  Editing line 1234...
🔨 Running test...
✅ Test passed!
```

## Agent Session Limits

**Important**: Each agent has API rate limits:
- If you see "Session limit reached", wait ~5 minutes
- Or launch a different agent type
- Or continue the work yourself

## How to Continue Agent Work

If an agent hits a rate limit:

```
Continue the work the agent was doing on [task].
Agent got to: [last step they completed]
Next steps: [what still needs to be done]
```

## Project-Specific Agent Commands

### Testing
```
Agent: Run full test suite and fix any failures
Command: python test_runner.py
Fix any issues found
Report: List of tests + pass/fail status
```

### Documentation
```
Agent: Document the [component]
Create: Implementation guide + usage examples
Format: Markdown with code samples
Save: docs/[component]_GUIDE.md
```

### Code Review
```
Agent: Review [file] for:
- Security issues
- Performance problems
- Code quality
- Best practices
Report findings with suggested fixes
```

## Tips for Maximum Productivity

### 1. Batch Similar Tasks
Instead of: "Fix bug A", then "Fix bug B", then "Fix bug C"
Do: "Launch 3 agents, each fix one bug simultaneously"

### 2. Use Context Files
Agents read:
- `.claude.md` - Quick project overview
- `CLAUDE_CONTEXT.md` - Detailed technical context
- `Todolist.md` - Full requirements

Reference these in agent prompts:
"Read .claude.md for context, then [task]"

### 3. Specify Output Format
"Create a report in markdown format"
"Generate test cases as Python code"
"Write implementation plan as numbered steps"

### 4. Set Success Criteria
"Implementation is complete when:
- All tests pass
- Feature works in extension
- Documentation is updated"

## Current Todo Agent Mapping

Based on your Todolist.md, here are suggested agent assignments:

### Ready to Launch Now:

**Agent 1: Image Size Controls** ⏳ Pending
- Complexity: Medium
- Files: 3 (content.js, popup.html, popup.js)
- Estimated time: 1 agent session

**Agent 2: 10 UI Variations** ⏳ Pending
- Complexity: Low (mostly CSS/HTML)
- Files: 10 new files
- Estimated time: 1 agent session

**Agent 3: Cache Deduplication** ⏳ Pending
- Complexity: Medium
- Files: server.js, mol2chemfig_server.py
- Estimated time: 1 agent session

**Agent 4: OPSIN 3D Integration** ⏳ Pending
- Complexity: High
- Needs: Research OPSIN API first
- Estimated time: 2 agent sessions

**Agent 5: New Compound Fix** ⏳ Pending
- Complexity: Medium
- Files: MoleculeViewer/server.js, content.js
- Estimated time: 1 agent session

### Already In Progress:

**PubChem Server** 🔄 In Progress
- Agent working on it
- May need continuation after rate limit

## Example Launch Command

To launch all pending agents in parallel:

```
Launch 5 agents in parallel to work on Todolist.md items:

1. general-purpose agent: Image size controls (todo #4)
2. general-purpose agent: 10 UI variations (todo #5)
3. general-purpose agent: Cache deduplication (todo #10)
4. general-purpose agent: OPSIN 3D integration (todo #9)
5. general-purpose agent: New compound fix (todo #11)

Each agent should:
- Read .claude.md for context
- Read Todolist.md for requirements
- Implement the feature completely
- Create tests
- Update documentation
- Report completion status
```

## Troubleshooting

### "Agent did something wrong"
Agents make mistakes! You can:
- Review their changes in git
- Ask me to fix what they did
- Revert and try again with clearer instructions

### "Agent got stuck"
If agent stops making progress:
- Interrupt and give more specific guidance
- Break task into smaller pieces
- Try a different agent type

### "Too many agents failed"
If multiple agents hit rate limits:
- Work on tasks yourself for 30 minutes
- Rate limits reset automatically
- Launch agents again later

## Next Steps

1. ✅ `.claude/settings.json` created - full auto-approval enabled
2. ✅ `.claude.md` created - quick reference
3. ✅ `CLAUDE_CONTEXT.md` updated - detailed context
4. ✅ TodoWrite updated - full detailed todos
5. ⏳ Ready to launch agents!

You can now say:
```
Launch agents to work on all pending todos in parallel.
Each agent should work autonomously and report results.
```

And I'll launch multiple agents to tackle your todolist simultaneously!
