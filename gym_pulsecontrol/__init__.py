from gym.envs.registration import register

register(id='PulseControl-v0',
    entry_point='gym_pulsecontrol.envs:PulseEnv',
)