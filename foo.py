from dataclasses import dataclass

first_bounce = True

@dataclass
class Electron(object):
    start_x: float
    start_y: float
    energy: float
    angle: float
    bounce_id: int

    initiating_electron: Electron
    secondaries: List

    def parent(self):
        return self.initiating_electron

    def gen_secondaries(self):


x = Electron(1,1,1,1,0,None)
y = Electron(1,1,1,1,0,x)


def generate_secondary(end_point, electron):
    # generates electrons with new start, energy... and
    # with bounce_id the same as electron.bounce_id + 1

    return secondaries


def track_electron_and_generate_secondarys(electron):
    end_point = track_electron(electron)
    secondaries = generate_secondary(end_point, electron)
    electron.secondaries = secondaries


electrons = []
if first_bounce:
    first_electron = gen_electron()
    electrons.append(first_electron)

all_electrons = []
while len(electrons) > 0:
    e = electrons.pop()
    new_e, bounce_id = track_electron_and_generate_secondarys(e)

    total_electron_count += len(new_e)

    all_electrons.extend(new_e)
    electrons.extend(new_e)


