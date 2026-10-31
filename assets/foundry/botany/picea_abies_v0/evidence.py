"""The evidence ledger, and the gate every generated number passes through.

Reading a parameter is the only way a botanical value reaches the generator, and
each read declares what it expects the row to be. A gap raises rather than
returning a zero, a superseded choice raises rather than quietly still applying,
and a value the fixture never declared it would consume raises rather than
becoming an untracked constant.
"""
from hashlib import sha256
from pathlib import Path
import json

MEASURED = ('verified-measurement', 'peer-reported-source-read')
CHOICE = 'explicit-calibration-choice'
QUALITATIVE = 'qualitative-observation'
FIXTURE_CLASS = 'visual-calibration-fixture'


class UnavailableMeasurement(LookupError):
    """A declared evidence gap was read. It carries no value and never will silently."""


class EvidenceError(ValueError):
    """The ledger or a read of it is malformed."""


class Ledger:
    def __init__(self, data, path=None):
        self.data = data
        self.path = path
        self.rows = {}
        for row in data['parameters']:
            if row['id'] in self.rows:
                raise EvidenceError(f'duplicate parameter id {row["id"]}')
            self.rows[row['id']] = row
        self._audit()
        self.sha256 = sha256(json.dumps(data, sort_keys=True, ensure_ascii=False).encode()).hexdigest()

    @classmethod
    def load(cls, path):
        path = Path(path)
        return cls(json.loads(path.read_text()), path)

    def _audit(self):
        if self.data.get('artifact_class') != FIXTURE_CLASS:
            raise EvidenceError('ledger is not marked a visual calibration fixture')
        if self.data.get('world_admission') is not False:
            raise EvidenceError('ledger must state world_admission false')
        classes = set(self.data['evidence_classes'])
        for pid, row in self.rows.items():
            if row['evidence_class'] not in classes:
                raise EvidenceError(f'{pid}: unknown evidence class {row["evidence_class"]}')
            gap = row['evidence_class'] == 'unavailable' or row.get('available') is False
            if gap and 'value' in row and row['value'] is not None:
                raise EvidenceError(f'{pid}: an unavailable measurement carries a value')
            if gap and not row.get('extraction_question'):
                raise EvidenceError(f'{pid}: a gap must name what would fill it')
            if row['evidence_class'] in MEASURED and not row.get('source'):
                raise EvidenceError(f'{pid}: a measurement must name its source')
            if row['evidence_class'] in MEASURED and row.get('source') not in self.data['sources']:
                raise EvidenceError(f'{pid}: source {row.get("source")} is not in the source list')
            if row['evidence_class'] == CHOICE and not row.get('rationale'):
                raise EvidenceError(f'{pid}: a calibration choice must say why')

    def row(self, pid):
        try:
            return self.rows[pid]
        except KeyError:
            raise EvidenceError(f'{pid} is not in the ledger') from None

    def value(self, pid, expect=None):
        row = self.row(pid)
        if row['evidence_class'] == 'unavailable' or row.get('available') is False:
            raise UnavailableMeasurement(f'{pid} is a declared evidence gap. {row["extraction_question"]}')
        if 'superseded_by' in row:
            raise EvidenceError(f'{pid} was superseded by {row["superseded_by"]}; read that instead')
        if expect and row['evidence_class'] not in expect:
            raise EvidenceError(f'{pid} is {row["evidence_class"]}, and this read expected {expect}')
        if 'value' not in row:
            raise EvidenceError(f'{pid} carries no value')
        return row['value']

    def measured(self, pid):
        return self.value(pid, MEASURED)

    def choice(self, pid):
        return self.value(pid, (CHOICE,))

    def classes(self):
        counts = {}
        for row in self.rows.values():
            counts[row['evidence_class']] = counts.get(row['evidence_class'], 0) + 1
        return counts

    def gaps(self):
        return {pid: row['extraction_question'] for pid, row in self.rows.items()
                if row['evidence_class'] == 'unavailable' or row.get('available') is False}


class ScopedLedger:
    """A ledger narrowed to what one fixture declared it would consume.

    The declaration lives in the fixture file, so a reviewer can read what the
    generator is allowed to touch without reading the generator.
    """

    def __init__(self, ledger, consumes, reports=()):
        self.ledger = ledger
        self.consumes = tuple(consumes)
        self.reports = tuple(reports)
        self.read = []
        unknown = [p for p in self.consumes + self.reports if p not in ledger.rows]
        if unknown:
            raise EvidenceError(f'fixture declares parameters absent from the ledger: {unknown}')

    def _gate(self, pid, allowed, what):
        if pid not in allowed:
            raise EvidenceError(f'{pid} is not in the fixture {what} list; declare it before reading it')
        if pid not in self.read:
            self.read.append(pid)

    def measured(self, pid):
        self._gate(pid, self.consumes, 'consumes')
        return self.ledger.measured(pid)

    def choice(self, pid):
        self._gate(pid, self.consumes, 'consumes')
        return self.ledger.choice(pid)

    def qualitative(self, pid):
        """Acknowledge a row that constrains a decision without supplying a number."""
        self._gate(pid, self.consumes, 'consumes')
        row = self.ledger.row(pid)
        if row['evidence_class'] not in (QUALITATIVE,) + MEASURED:
            raise EvidenceError(f'{pid} is not an observation')
        return row.get('value')

    def report(self, pid):
        self._gate(pid, self.consumes + self.reports, 'consumes or reports')
        return self.ledger.value(pid)

    def unread(self):
        return [p for p in self.consumes if p not in self.read]
