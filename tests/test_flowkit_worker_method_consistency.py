from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FLOWKIT_PYTHON = ROOT / ".tool_envs" / "flowkit" / "Scripts" / "python.exe"


@pytest.mark.skipif(not FLOWKIT_PYTHON.is_file(), reason="locked FlowKit environment unavailable")
def test_flowjo_method_consistency_includes_membership_logic_without_fcs_inputs():
    script = r'''
import importlib.util
import sys
import tempfile
from pathlib import Path

root = Path(sys.argv[1])
sys.path.insert(0, str(root / "integrations" / "workers"))
spec = importlib.util.spec_from_file_location(
    "flowkit_worker", root / "integrations" / "workers" / "flowkit_worker.py"
)
worker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(worker)

class Dimension:
    id = "FL1-A"

class Sample:
    def __init__(self, sample_id):
        self.sample_id = sample_id
    def get_metadata(self):
        return {"fil": f"{self.sample_id}.fcs"}

class Workspace:
    def __init__(self, gates):
        self.samples = {key: Sample(key) for key in gates}
        self.gates = gates
    def get_sample_ids(self):
        return list(self.samples)
    def get_sample(self, sample_id):
        return self.samples[sample_id]
    def get_gate_ids(self, sample_id):
        return [("Live", ("root",))]
    def get_gate(self, sample_id, gate_name, gate_path):
        return self.gates[sample_id]
    def get_transforms(self, sample_id):
        return {}
    def get_comp_matrix(self, sample_id):
        return None

def rectangle(complement=False):
    gate_type = type("RectangleGate", (), {})
    gate = gate_type()
    gate.dimensions = [Dimension()]
    gate.use_complement = complement
    return gate

def boolean(operation, reference):
    gate_type = type("BooleanGate", (), {})
    gate = gate_type()
    gate.dimensions = None
    gate.type = operation
    gate.gate_refs = [{"ref": reference, "path": ("root",), "complement": False}]
    return gate

def compare(gates):
    cached = [
        {"keywords": {"fil": f"{sample_id}.fcs"}, "counts": {
            (("root",), "Live"): {"source_gate_type": type(gate).__name__}
        }}
        for sample_id, gate in gates.items()
    ]
    input_ids = {sample_id: f"{sample_id}.fcs" for sample_id in gates}
    return worker._verify_wsp_gate_method(Workspace(gates), cached, input_ids)[0]

assert compare({"H1": rectangle(), "H2": rectangle()})["status"] == "CONSISTENT"
assert compare({"H1": rectangle(), "H2": rectangle(True)})["status"] == "INCONSISTENT"
assert compare({"H1": boolean("and", "Live"), "H2": boolean("or", "Live")})["status"] == "INCONSISTENT"
assert compare({"H1": boolean("and", "Live"), "H2": boolean("and", "Other")})["status"] == "INCONSISTENT"

wsp_xml = """<Workspace xmlns:g="http://www.isac-net.org/std/Gating-ML/v2.0/gating"
    xmlns:d="http://www.isac-net.org/std/Gating-ML/v2.0/datatypes">
  <Sample><SampleNode sampleID="1" name="H1.fcs"><Subpopulations>
    <Population name="Live"><Gate><g:EllipsoidGate>
      <g:dimension g:transformation-ref="biex">
        <d:fcs-dimension d:name="FL1-A" />
      </g:dimension>
      <g:dimension><d:fcs-dimension d:name="SSC-A" /></g:dimension>
    </g:EllipsoidGate></Gate></Population>
  </Subpopulations></SampleNode><Transformations /></Sample>
</Workspace>"""
transform_xml = """<Workspace xmlns:t="http://www.isac-net.org/std/Gating-ML/v2.0/transformations"
    xmlns:d="http://www.isac-net.org/std/Gating-ML/v2.0/datatypes">
  <Sample><SampleNode sampleID="1" name="H1.fcs" />
    <Transformations>
      <hyperlog><d:parameter d:name="FL1-A" /></hyperlog>
      <t:biex t:neg="0.1" t:width="0.5" t:length="512" t:maxRange="262144" t:pos="4.5">
        <d:parameter d:name="RL1-A" />
      </t:biex>
      <t:logicle t:T="262144" t:W="0.5" t:M="4.5">
        <d:parameter d:name="VL1-A" />
      </t:logicle>
      <t:fasinh t:T="262144" t:A="0" t:M="NaN">
        <d:parameter d:name="WL1-A" />
      </t:fasinh>
    </Transformations>
  </Sample>
</Workspace>"""
valid_transform_xml = """<Workspace xmlns:t="http://www.isac-net.org/std/Gating-ML/v2.0/transformations"
    xmlns:d="http://www.isac-net.org/std/Gating-ML/v2.0/datatypes">
  <Sample><SampleNode sampleID="1" name="H1.fcs" />
    <Transformations>
      <t:logicle t:T="262144" t:W="0.5" t:M="4.5" t:A="0">
        <d:parameter d:name="FL1-A" />
      </t:logicle>
    </Transformations>
  </Sample>
</Workspace>"""
with tempfile.TemporaryDirectory() as temp_dir:
    wsp_path = Path(temp_dir) / "synthetic.wsp"
    wsp_path.write_text(wsp_xml, encoding="utf-8")
    cached = [{"workspace_sample_id": "__wsp_sample_1_0"}]
    nodes = worker._wsp_transformed_ellipsoid_nodes(wsp_path, cached)
    assert len(nodes) == 1
    assert nodes[0]["reason_code"] == "FLOWJO_ELLIPSOID_TRANSFORM_UNSUPPORTED"
    assert nodes[0]["channels"] == ["FL1-A"]
    assert nodes[0]["parent_path"] == ["root"]
    assert worker._wsp_unsupported_gate_nodes(wsp_path, cached) == []
    result = worker._flowjo_unsupported_gate_result(
        wsp_path, cached, {"__wsp_sample_1_0": "H1.fcs"}
    )
    assert result["status"] == "WAITING_FOR_USER"
    assert result["unresolved_nodes"][0]["sample_id"] == "H1.fcs"

    wsp_path.write_text(wsp_xml.replace(' g:transformation-ref="biex"', ""), encoding="utf-8")
    assert worker._wsp_transformed_ellipsoid_nodes(wsp_path, cached) == []
    assert worker._flowjo_unsupported_gate_result(wsp_path, cached, {}) is None

    unsupported_gate_xml = """<Workspace xmlns:g="http://www.isac-net.org/std/Gating-ML/v2.0/gating"
        xmlns:d="http://www.isac-net.org/std/Gating-ML/v2.0/datatypes">
      <Sample><SampleNode sampleID="1" name="H1.fcs"><Subpopulations>
        <Population name="Live"><Gate><g:RectangleGate /></Gate><Subpopulations>
          <Population name="Quadrants"><Gate><g:QuadrantGate>
            <g:dimension><d:fcs-dimension d:name="FL1-A" /></g:dimension>
            <g:dimension><d:fcs-dimension d:name="SSC-A" /></g:dimension>
          </g:QuadrantGate></Gate></Population>
        </Subpopulations></Population>
        <Population name="Custom"><Gate><g:SectorGate /></Gate></Population>
      </Subpopulations></SampleNode><Transformations /></Sample>
    </Workspace>"""
    wsp_path.write_text(unsupported_gate_xml, encoding="utf-8")
    gate_nodes = worker._wsp_unsupported_gate_nodes(wsp_path, cached)
    assert {node["reason_code"] for node in gate_nodes} == {"FLOWJO_GATE_TYPE_UNSUPPORTED"}
    assert {node["source_gate_type"] for node in gate_nodes} == {"QuadrantGate", "SectorGate"}
    quadrant_node = next(node for node in gate_nodes if node["name"] == "Quadrants")
    assert quadrant_node["channels"] == ["FL1-A", "SSC-A"]
    assert quadrant_node["parent_path"] == ["root", "Live"]
    gate_result = worker._flowjo_unsupported_gate_result(
        wsp_path, cached, {"__wsp_sample_1_0": "H1.fcs"}
    )
    assert gate_result["status"] == "WAITING_FOR_USER"
    assert all(node["sample_id"] == "H1.fcs" for node in gate_result["unresolved_nodes"])

    wsp_path.write_text(transform_xml, encoding="utf-8")
    transform_nodes = worker._wsp_transform_nodes(wsp_path, cached)
    assert {node["reason_code"] for node in transform_nodes} == {
        "FLOWJO_TRANSFORM_UNSUPPORTED", "FLOWJO_TRANSFORM_PARAMETER_UNSUPPORTED",
        "FLOWJO_TRANSFORM_PARAMETERS_INVALID",
    }, transform_nodes
    assert {node["transform_type"] for node in transform_nodes} == {
        "hyperlog", "biex", "logicle", "fasinh"
    }, transform_nodes
    transform_result = worker._flowjo_unsupported_gate_result(
        wsp_path, cached, {"__wsp_sample_1_0": "H1.fcs"}
    )
    assert transform_result["status"] == "WAITING_FOR_USER"
    assert all(node["sample_id"] == "H1.fcs" for node in transform_result["unresolved_nodes"])

    wsp_path.write_text(valid_transform_xml, encoding="utf-8")
    assert worker._wsp_transform_nodes(wsp_path, cached) == []

    missing_transforms_xml = """<Workspace><Sample>
      <SampleNode sampleID="1" name="H1.fcs" />
    </Sample></Workspace>"""
    wsp_path.write_text(missing_transforms_xml, encoding="utf-8")
    missing_nodes = worker._wsp_transform_nodes(wsp_path, cached)
    assert missing_nodes[0]["reason_code"] == "FLOWJO_TRANSFORM_BLOCK_MISSING", missing_nodes
'''
    completed = subprocess.run(
        [str(FLOWKIT_PYTHON), "-c", script, str(ROOT)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
