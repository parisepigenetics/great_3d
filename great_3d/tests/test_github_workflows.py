"""
Unit tests for GitHub Actions workflow configurations.
"""

import os
import pytest

# Check if yaml module is available
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

yaml = pytest.importorskip("yaml")


@pytest.fixture
def opencode_workflow():
    """Load the opencode workflow YAML file."""
    workflow_path = os.path.join(
        os.path.dirname(__file__),
        '..',
        '..',
        '.github',
        'workflows',
        'opencode.yml'
    )
    with open(workflow_path, 'r') as f:
        return yaml.safe_load(f)


def test_opencode_triggers_on_issue_comment(opencode_workflow):
    """Test that opencode workflow triggers on issue_comment with created type."""
    assert 'on' in opencode_workflow
    assert 'issue_comment' in opencode_workflow['on']
    assert opencode_workflow['on']['issue_comment']['types'] == ['created']


def test_opencode_triggers_on_pull_request_review_comment(opencode_workflow):
    """Test that opencode workflow triggers on pull_request_review_comment with created type."""
    assert 'on' in opencode_workflow
    assert 'pull_request_review_comment' in opencode_workflow['on']
    assert opencode_workflow['on']['pull_request_review_comment']['types'] == ['created']


def test_opencode_triggers_with_specific_commands(opencode_workflow):
    """Test that opencode workflow only runs when specific commands are present in comments."""
    opencode_job = opencode_workflow['jobs']['opencode']
    
    assert 'if' in opencode_job
    
    # The conditional should check for /oc or /opencode commands
    condition = opencode_job['if']
    assert '/oc' in condition
    assert '/opencode' in condition
    assert 'contains(github.event.comment.body' in condition
    assert 'startsWith(github.event.comment.body' in condition


def test_opencode_runs_on_ubuntu_latest(opencode_workflow):
    """Test that opencode workflow runs on ubuntu-latest."""
    opencode_job = opencode_workflow['jobs']['opencode']
    assert opencode_job['runs-on'] == 'ubuntu-latest'


def test_opencode_has_correct_permissions(opencode_workflow):
    """Test that opencode workflow has correct permissions configured."""
    opencode_job = opencode_workflow['jobs']['opencode']
    
    assert 'permissions' in opencode_job
    permissions = opencode_job['permissions']
    
    # Verify all required permissions
    assert permissions['id-token'] == 'write'
    assert permissions['contents'] == 'read'
    assert permissions['pull-requests'] == 'read'
    assert permissions['issues'] == 'read'


def test_opencode_uses_anomalyco_action(opencode_workflow):
    """Test that opencode workflow uses anomalyco/opencode/github action."""
    opencode_job = opencode_workflow['jobs']['opencode']
    
    # Find the step that runs opencode
    opencode_step = None
    for step in opencode_job['steps']:
        if step.get('name') == 'Run opencode':
            opencode_step = step
            break
    
    assert opencode_step is not None, "Run opencode step not found"
    assert opencode_step['uses'] == 'anomalyco/opencode/github@latest'


def test_opencode_provides_api_key_environment_variable(opencode_workflow):
    """Test that opencode workflow provides OPENCODE_API_KEY environment variable."""
    opencode_job = opencode_workflow['jobs']['opencode']
    
    # Find the step that runs opencode
    opencode_step = None
    for step in opencode_job['steps']:
        if step.get('name') == 'Run opencode':
            opencode_step = step
            break
    
    assert opencode_step is not None, "Run opencode step not found"
    assert 'env' in opencode_step
    assert 'OPENCODE_API_KEY' in opencode_step['env']
    
    # Verify it references the secret
    api_key = opencode_step['env']['OPENCODE_API_KEY']
    assert 'secrets.OPENCODE_API_KEY' in api_key
