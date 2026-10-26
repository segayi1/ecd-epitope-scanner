from html.parser import HTMLParser
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ResourceParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.resources = []

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name in {"src", "href"} and value:
                self.resources.append((tag, name, value))


def test_offline_template_layout_and_interaction_contract():
    template = (ROOT / "web/offline_template.html").read_text(encoding="utf-8")

    assert "--bg:#f3f6f7" in template
    assert "--surface:#fff" in template
    assert 'grid-template-columns:minmax(0,2fr) minmax(340px,1fr)' in template
    assert '#main[data-step="5"],#main[data-step="6"],#main[data-step="7"]' in template
    assert 'var STEPS = ["step.surface","step.tolerance","step.regions","step.interface","step.design","step.custom","step.overview"]' in template
    assert 'aria-current=' in template and 'data-step=' in template
    assert 'sequenceExpanded: false' in template
    assert 'aria-controls="seq-content"' in template and 'sequenceExpanded = !state.sequenceExpanded' in template
    assert 'id="btn-legend"' in template and 'aria-expanded="false"' in template
    assert 'id="cp-wt-aa"' in template and 'id="cp-mu-aa"' in template
    assert 'state.lastPick === pos' in template and 'applyStyles(); renderSeq();' in template
    assert 'createRecommendationRenderer' in template
    design_step = template.split('function stepDesign(){', 1)[1].split('\nfunction ', 1)[0]
    assert 'var designNote = recommendations.length' in design_step
    assert '各分区按设计原则独立给出通过门槛的推荐，数量不受全局精简面板预算限制。' in design_step
    assert 'Each region independently recommends combinations that pass the gates' in design_step
    assert ': t("design.note")' in design_step

    parser = ResourceParser()
    parser.feed(template)
    assert parser.resources == []

