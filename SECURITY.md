# Security Policy

## Reporting a vulnerability

If you discover a security issue (for example in deployment scripts, a web
interface, or CI workflows), please email aaron.weiskittel@maine.edu rather than
opening a public issue.

## Scope

This repository contains scientific research code, not a production service.
Where the repository includes deployment or pipeline infrastructure that may run
in a server context, security reports for those components are welcome.

## Credentials

Never commit credentials, API tokens, SSH keys, or passwords to this repository.
Check `.gitignore` before committing, and if a secret is ever exposed, revoke and
rotate it immediately.
