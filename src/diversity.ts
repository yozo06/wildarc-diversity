/**
 * Biodiversity indices from species abundances.
 *
 * Sources:
 *   Magurran, A. E. (2004). Measuring Biological Diversity. Blackwell Publishing. (ch. 4)
 *   Jost, L. (2006). Entropy and diversity. Oikos 113(2):363–375.
 *   Hurlbert, S. H. (1971). The nonconcept of species diversity: a critique and alternative parameters.
 *     Ecology 52(4):577–586.
 */

/** A computed value, or an explicit refusal when the method does not apply. Never NaN, Infinity or -0. */
export type Measure = { value: number } | { value: null; reason: string };

/** Abundance per species: an array of counts, or a record keyed by species name. */
export type Abundances = readonly number[] | Readonly<Record<string, number>>;

export interface DiversityOptions {
    /** Logarithm base for `shannon` only (default Math.E, i.e. nats). No other field depends on it. */
    logBase?: number;
}

export interface DiversityProfile {
    /** N: sum of the abundances greater than 0. */
    totalAbundance: number;
    /** S: number of species with an abundance greater than 0. */
    richness: number;
    /** Shannon index H′ = −Σ pᵢ log_b pᵢ. */
    shannon: Measure;
    /** Pielou's evenness J = H′ / ln S. Does not depend on the log base. */
    pielouEvenness: Measure;
    /** Simpson's λ = Σ nᵢ(nᵢ−1) / (N(N−1)), the unbiased finite-sample form. */
    simpsonLambda: Measure;
    /** 1 − λ with the unbiased λ: the chance that two individuals drawn without replacement differ in species. */
    giniSimpson: Measure;
    /** Hill numbers (effective number of species): q0 = S, q1 = exp(H′), q2 = 1/Σpᵢ². */
    hill: { q0: number; q1: Measure; q2: Measure };
}

const EMPTY_COMMUNITY = 'empty community: no individuals counted';
const TOO_FEW_SPECIES = 'evenness is undefined for fewer than 2 species';
const TOO_FEW_INDIVIDUALS = "Simpson's index needs at least 2 individuals";
const NOT_WHOLE_COUNTS = "Simpson's unbiased index needs whole-number counts";

// -0 is normalised to 0: JSON.stringify hides the sign, Object.is and 1/x do not.
const computed = (value: number): Measure => ({ value: value === 0 ? 0 : value });
const refused = (reason: string): Measure => ({ value: null, reason });

function toEntries(abundances: Abundances): [string, number][] {
    if (Array.isArray(abundances)) {
        return Array.from(abundances as readonly number[], (n, i) => [`index ${i}`, n]);
    }
    return Object.entries(abundances as Readonly<Record<string, number>>);
}

export function diversityProfile(abundances: Abundances, options: DiversityOptions = {}): DiversityProfile {
    const logBase = options.logBase ?? Math.E;
    if (!Number.isFinite(logBase) || logBase <= 0 || logBase === 1) {
        throw new RangeError(`logBase must be a finite number greater than 0 and not 1, got ${logBase}`);
    }

    // Species with abundance 0 are absent: they count towards neither S nor N.
    const counts: number[] = [];
    for (const [label, n] of toEntries(abundances)) {
        if (!Number.isFinite(n) || n < 0) {
            throw new RangeError(`abundance for ${label} must be a finite number >= 0, got ${n}`);
        }
        if (n > 0) counts.push(n);
    }
    const totalAbundance = counts.reduce((sum, n) => sum + n, 0);
    const richness = counts.length;
    if (!Number.isFinite(totalAbundance)) {
        throw new RangeError('total abundance is too large to represent as a finite number');
    }

    if (richness === 0) {
        return {
            totalAbundance: 0,
            richness: 0,
            shannon: refused(EMPTY_COMMUNITY),
            pielouEvenness: refused(EMPTY_COMMUNITY),
            simpsonLambda: refused(EMPTY_COMMUNITY),
            giniSimpson: refused(EMPTY_COMMUNITY),
            hill: { q0: 0, q1: refused(EMPTY_COMMUNITY), q2: refused(EMPTY_COMMUNITY) },
        };
    }

    // Shannon index H′ = −Σ pᵢ ln pᵢ, in natural logs (Magurran 2004, ch. 4).
    // Σ pᵢ² is kept for the Hill number of order 2.
    let shannonNats = 0;
    let sumPSquared = 0;
    for (const n of counts) {
        const p = n / totalAbundance;
        // p ln p → 0 as p → 0, so a share too small to represent adds nothing (0 × ln 0 would be NaN).
        if (p > 0) shannonNats -= p * Math.log(p);
        sumPSquared += p * p;
    }

    // Pielou's evenness J = H′ / ln S (Magurran 2004, ch. 4). ln 1 = 0, so it is undefined for one species.
    const pielouEvenness = richness < 2 ? refused(TOO_FEW_SPECIES) : computed(shannonNats / Math.log(richness));

    // Simpson's λ = Σ nᵢ(nᵢ−1) / (N(N−1)): the chance that two individuals drawn without replacement
    // are the same species (Magurran 2004, ch. 4). Only meaningful for counts of individuals.
    let simpsonLambda: Measure;
    if (!counts.every((n) => Number.isInteger(n))) {
        simpsonLambda = refused(NOT_WHOLE_COUNTS);
    } else if (totalAbundance < 2) {
        simpsonLambda = refused(TOO_FEW_INDIVIDUALS);
    } else {
        // Summed as (nᵢ/N)·((nᵢ−1)/(N−1)) so that very large counts cannot overflow.
        const lambda = counts.reduce((sum, n) => sum + (n / totalAbundance) * ((n - 1) / (totalAbundance - 1)), 0);
        simpsonLambda = computed(lambda);
    }

    // 1 − λ: the chance that two individuals drawn without replacement are different species. With the
    // unbiased λ this is (N/(N−1))(1 − Σpᵢ²), Hurlbert's (1971) probability of interspecific encounter. The
    // Gini-Simpson index in Jost (2006, Table 1) is the plug-in form 1 − Σpᵢ², which is 1 − 1/hill.q2.
    const giniSimpson = simpsonLambda.value === null ? refused(simpsonLambda.reason) : computed(1 - simpsonLambda.value);

    return {
        totalAbundance,
        richness,
        shannon: computed(shannonNats / Math.log(logBase)),
        pielouEvenness,
        simpsonLambda,
        giniSimpson,
        // Hill numbers, the "effective number of species" (Jost 2006): q0 = S, q1 = exp(H′), q2 = 1/Σpᵢ².
        hill: { q0: richness, q1: computed(Math.exp(shannonNats)), q2: computed(1 / sumPSquared) },
    };
}
