import cobra
import cobra.test
import gym
from gym import spaces
from gym.utils import seeding
import numpy as np
import numpy.random as npr
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import copy
import utils
from joblib import Parallel

class FBA_Env(gym.Env):      
    def __init__(self, model, df):
        self.model = model.copy()
        self.cur_model = self.model.copy()
        self.rxns = copy.deepcopy(model.reactions)
        self.n_rxns = len(self.rxns)
        # TODO Add metabs
        self.metabs = copy.deepcopy(model.metabolites)
        n_metabs = len(self.metabs)
        # LEARNING FROM EXPERIMENTAL DATA:
        self.df = df
        self.cond_scores = df['AVG.1']
        # Add 1 for not removing any
        self.action_space = spaces.Discrete(self.n_rxns) # TODO: Think about using MultiBinary
        self.action_types = ('remove', 'add')
        #spaces.Dict({"reaction": spaces.Discrete(n_rxns + 1), "metabolite": spaces.Discrete(n_metabs + 1)})
        # Whether or not a reaction is present
        self.observation_space = spaces.MultiBinary(self.n_rxns)
        # In case they forget to reset
        self.state = np.ones(self.n_rxns, dtype=np.int8)
        self.time = 0
        
        self.unique_fluxes = {}
        
        self.log = []
        
        self._seed()
        #spaces.Dict({"reaction": spaces.Discrete(2), "metabolite": spaces.Discrete(3)})
        
    def _seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]

    def _get_obs(self):
        return self.state
    
    def _take_action(self, state, action):
        # Returns new state
        return NotImplemented
    
    def _evaluate(self, state_rxns):
        return NotImplemented
    
    def reset(self):
        # TODO, choose random starting state
        self.state = np.ones(self.n_rxns, dtype=np.int8)
        self.cur_model = self.model.copy()
        self.last_action = None
        self.time = 0
        
        return self._get_obs()

    def step(self,action):
        # TODO: do better than this?
        self._take_action(self.state, action)
        #self.state = new_state
        reward = self._evaluate(self.state)

        done = False
        #if self.time > 300:
        #    done = True

        self.time += 1
            
        return self.rxns[action[0]], reward, done, action

# TODO: Look into GoalEnv
class FBA_Step_Env(FBA_Env):  
    def _take_action(self, state, act):
        action, action_type = act
        if action_type == 'add':
            self.state[action] = 1
            rxn = self.rxns[action]
            self.cur_model.add_reaction(rxn)
        else:
        # Returns new state
            self.state[action] = 0
            rxn = self.rxns[action]
            if rxn in self.cur_model.reactions:
                self.cur_model.reactions.get_by_id(rxn.id).remove_from_model()
            else:
                self.log.append(('rxn {0} not in model with state {1}'.format(rxn, self.state)))
        #return state
    
    def _evaluate_slim(self, state_rxns):
        print 'N RXNS: {0}'.format(len(self.cur_model.reactions))
        objs = utils.add_addl_reactants(self.cur_model, self.df)
        if not tuple(objs) in self.unique_fluxes:
            self.unique_fluxes[tuple(objs)] = self.state.copy()
        if sum(objs) < 0.01 or np.isnan(objs).any() or not objs:
            return -1000
        corr = scipy.stats.spearmanr(objs, self.cond_scores)
        return corr
    
    def _evaluate(self, state_rxns):
        #objs = utils.add_addl_reactants(self.cur_model, self.df)
        objs, fluxes = utils.gen_fluxes_addl_reactants(self.cur_model, self.df, self.parallel)
        if not tuple(objs) in self.unique_fluxes:
            self.unique_fluxes[tuple(objs)] = (fluxes, state_rxns)
        if sum(objs) < 0.01 or np.isnan(objs).any() or not objs:
            return -1000
        corr = scipy.stats.spearmanr(objs, self.cond_scores)
        return corr
    
    def _evaluate_flux(self, state_rxns):
        #objs = utils.add_addl_reactants(self.cur_model, self.df)
        obj, fluxes = utils.gen_fluxes(self.cur_model)
        if not tuple(fluxes) in self.unique_fluxes:
            self.unique_fluxes[tuple(fluxes)] = state_rxns
        if obj < 10 ** -7:
            return -1000
        return obj

class RandomAgent(object):
    def __init__(self, action_space):
        self.action_space = action_space

    def act(self, observation, reward, done, prev_action):
        return (self.action_space.sample(), 'remove')

class AddBackAgent(object):
    def __init__(self, action_space):
        self.action_space = action_space

    def act(self, observation, reward, done, prev_action):
        if reward < -1:
            return (prev_action[0], 'add')
        return (self.action_space.sample(), 'remove')

if __name__ == '__main__':
    model = cobra.io.read_sbml_model('../models/ecoli_cf_base.sbml')

    df = pd.read_csv('../data/Karim_MetEng_2018_Figure2_Data.csv')
    df.drop(columns=['Area_1', 'Area_2', 'Conc_1', 'Conc_2'], inplace=True)
    df.head()

    env = FBA_Step_Env(model, df)
    agent = AddBackAgent(env.action_space)
    max_reward = ([0], None)
    min_reward = ([0], None)
    prev_action, done = None, False
    for i_episode in range(1):
        print 'Episode {0}'.format(i_episode)
        observation = env.reset()
        reward = 0
        for t in range(3):
            #env.render()
            #print(observation)
            action = agent.act(observation, reward, done, prev_action)
            observation, reward, done, prev_action = env.step(action)
            if np.isnan(reward).any():
                print 'NAN'
                break
            if t % 10 == 0:
                print 'Time {0}, reward: {1}'.format(t, reward)
            if reward[0] > max_reward[0][0]:
                max_reward = (reward, env.state)
            if reward[0] < min_reward[0][0]:
                min_reward = (reward, env.state)
            if done:
                print("Episode finished after {} timesteps".format(t+1))
                break
    print max_reward, min_reward

    dis_rxns = []
    obj_l = []
    i = 0
    for objs, rxns in env.unique_fluxes.items():
        i += 1
        if sum(objs) < 0.01:
            continue
        dis_rxns.append(rxns)
        obj_l.append(objs)
    print i
    print len(dis_rxns)

    max_len = (0, None)
    for rxn in range(len(dis_rxns)):
	for flux_ser in dis_rxns[rxn][0]:
	    if flux_ser.shape[0] > max_len[0]:
		max_len = (flux_ser.shape[0], flux_ser)
    ind = max_len[1].index

    fluxes = []
    for rxn in range(len(dis_rxns)):
        experiments = []
        for flux_ser in dis_rxns[rxn][0]:
            flux_ser_pad = flux_ser.reindex(ind)
            experiments.append(np.array(flux_ser_pad))
        experiment = np.stack(experiments, axis=0)
        fluxes.append(experiment)
    flux_arr = np.stack(fluxes, axis=0)
    flux_arr.shape


    np.save('../data/fluxes_ecoli_but_test', flux_arr)
    np.save('../data/rxns_ecoli_but_test', env.rxns)
