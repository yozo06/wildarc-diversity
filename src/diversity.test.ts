import { describe, it, expect } from 'vitest';
import { diversityProfile } from './index';
import type { Abundances, Measure } from './index';

const v = (m: Measure) => m.value;

describe('diversityProfile', () => {
    it('even community [10,10,10,10]', () => {
        const p = diversityProfile([10, 10, 10, 10]);
        expect(p.totalAbundance).toBe(40);
        expect(p.richness).toBe(4);
        expect(v(p.shannon)).toBeCloseTo(1.3862943611, 10); // ln 4
        expect(v(p.pielouEvenness)).toBeCloseTo(1, 10);
        expect(v(p.simpsonLambda)).toBeCloseTo(0.2307692308, 10); // 360/1560
        expect(v(p.giniSimpson)).toBeCloseTo(0.7692307692, 10);
        expect(p.hill.q0).toBe(4);
        expect(v(p.hill.q1)).toBeCloseTo(4, 10);
        expect(v(p.hill.q2)).toBeCloseTo(4, 10);
    });

    it('uneven community [50,30,20]', () => {
        const p = diversityProfile([50, 30, 20]);
        expect(v(p.shannon)).toBeCloseTo(1.0296530141, 10);
        expect(v(p.pielouEvenness)).toBeCloseTo(0.9372305632, 10);
        expect(v(p.simpsonLambda)).toBeCloseTo(0.3737373737, 10); // 3700/9900
        expect(v(p.giniSimpson)).toBeCloseTo(0.6262626263, 10);
        expect(v(p.hill.q1)).toBeCloseTo(2.8000940729, 10);
        expect(v(p.hill.q2)).toBeCloseTo(2.6315789474, 10); // 1/0.38
    });

    it('record input equals array input', () => {
        expect(diversityProfile({ 'Grevillea robusta': 50, 'Coffea arabica': 30, 'Mangifera indica': 20 }))
            .toEqual(diversityProfile([50, 30, 20]));
    });

    it('logBase changes shannon only', () => {
        const p = diversityProfile([10, 10, 10, 10], { logBase: 2 });
        expect(v(p.shannon)).toBeCloseTo(2, 10);
        expect(v(p.hill.q1)).toBeCloseTo(4, 10);
        expect(v(p.pielouEvenness)).toBeCloseTo(1, 10);
    });

    it('zero counts do not add species', () => {
        expect(diversityProfile([10, 0, 10, 0]).richness).toBe(2);
    });

    it('single species: shannon is +0, evenness null', () => {
        const p = diversityProfile([5]);
        expect(p.shannon.value).toBe(0); // fails on -0 (Object.is)
        expect(p.pielouEvenness).toEqual({ value: null, reason: 'evenness is undefined for fewer than 2 species' });
        expect(v(p.simpsonLambda)).toBe(1);
        expect(v(p.hill.q2)).toBe(1);
    });

    it('one individual: simpson null', () => {
        expect(diversityProfile([1]).simpsonLambda)
            .toEqual({ value: null, reason: "Simpson's index needs at least 2 individuals" });
    });

    it('two singletons [1,1]: lambda is exactly 0, gini 1', () => {
        const p = diversityProfile([1, 1]);
        expect(v(p.simpsonLambda)).toBe(0);
        expect(v(p.giniSimpson)).toBe(1);
    });

    it('non-integer abundances: shannon/hill computed, simpson null', () => {
        const p = diversityProfile([2.5, 1.5]);
        expect(v(p.shannon)).toBeCloseTo(0.6615632382, 10);
        expect(v(p.hill.q2)).toBeCloseTo(1.8823529412, 10);
        expect(p.simpsonLambda.value).toBeNull();
        expect(p.giniSimpson).toEqual({ value: null, reason: "Simpson's unbiased index needs whole-number counts" });
    });

    it.each([[[]], [[0, 0]], [{}]])('empty community %j: every measure null with reason', (input) => {
        const p = diversityProfile(input as Abundances);
        expect(p.totalAbundance).toBe(0);
        expect(p.richness).toBe(0);
        expect(p.hill.q0).toBe(0);
        for (const m of [p.shannon, p.pielouEvenness, p.simpsonLambda, p.giniSimpson, p.hill.q1, p.hill.q2]) {
            expect(m).toEqual({ value: null, reason: 'empty community: no individuals counted' });
        }
    });

    it.each([[[-1]], [[NaN]], [[Infinity]], [{ teak: -2 }]])('invalid abundance %j throws RangeError', (input) => {
        expect(() => diversityProfile(input as Abundances)).toThrow(RangeError);
    });

    it('RangeError names the offending entry', () => {
        expect(() => diversityProfile({ teak: -2 })).toThrow(/teak/);
        expect(() => diversityProfile([3, -1])).toThrow(/index 1/);
    });

    it.each([0, 1, -2, NaN, Infinity])('invalid logBase %s throws RangeError', (b) => {
        expect(() => diversityProfile([1, 2], { logBase: b })).toThrow(RangeError);
    });

    // Extreme magnitudes: the plain formulas overflow or underflow to NaN here.
    it('total abundance beyond the largest finite number throws RangeError', () => {
        expect(() => diversityProfile([1e308, 1e308])).toThrow(RangeError);
    });

    it('a species too rare to register (p underflows to 0) contributes 0, not NaN', () => {
        const p = diversityProfile([1e-300, 1e300]);
        expect(p.richness).toBe(2);
        expect(v(p.shannon)).toBe(0);
        expect(v(p.pielouEvenness)).toBe(0);
        expect(v(p.hill.q1)).toBe(1);
    });

    it('very large whole-number counts keep simpson finite', () => {
        const p = diversityProfile([1e200, 1e200]);
        expect(v(p.simpsonLambda)).toBeCloseTo(0.5, 10);
        expect(v(p.giniSimpson)).toBeCloseTo(0.5, 10);
    });
});
