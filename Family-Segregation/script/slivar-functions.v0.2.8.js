function x_hom(s) {
  if (!s) return false;
  if (s.sex == "male" || s.sex == 1) {
    return s.alts >= 1; // Male hemizygote
  }
  if (s.sex == "female" || s.sex == 2) {
    return s.hom_alt || s.alts == 2; // Female homozygote
  }

  return s.hom_alt;
}

function denovo(kid, mom, dad) {
  // Check that child is het (0/1) and both parents are hom_ref (0/0)
  return kid.het && mom.hom_ref && dad.hom_ref;
}

function x_denovo(kid, mom, dad) {
  if (kid.sex != 'male') { return false; }
  return kid.alts >= 1 && mom.hom_ref && dad.hom_ref;
}

function uniparent_disomy(kid, mom, dad) {
  return (kid.alts == 0 || kid.alts == 2) && 
         ((mom.alts == 2 && dad.alts == 0) || (mom.alts == 0 && dad.alts == 2));
}

function recessive(kid, mom, dad) {
  return kid.hom_alt && mom.het && dad.het;
}

function x_recessive(kid, mom, dad) { 
  return kid.sex == 'male' && mom.het && dad.hom_ref && kid.alts >= 1;
}

// Heterozygous (1 side of compound het)
function solo_ch_het_side(sample) {
  return sample.het;
}

function comphet_side(kid, mom, dad) {
  return kid.het && 
         (solo_ch_het_side(mom) != solo_ch_het_side(dad)) && 
         mom.alts != 2 && dad.alts != 2;
}

// Assume that mom and kid are affected
function fake_auto_dom(kid, mom, dad) {
  return kid.het && mom.het && dad.hom_ref;
}

function segregating_dominant_x(s) {
  if (!s) return false;

  // UNAFFECTED INDIVIDUALS (Parents, Siblings, Relatives)
  // Must NOT carry an alternate allele (allows 0/0 or missing ./. call)
  if (!s.affected) {
    return s.hom_ref || s.alts < 0;
  }

  // AFFECTED FEMALES
  if (s.sex == "female" || s.sex == 2 || s.sex == "2") {
    // In X-dominant, affected females MUST be heterozygous
    if (!s.het) return false;

    // Transmission Check (if parents are present & sequenced)
    var has_mom = ("mom" in s && s.mom && s.mom.alts >= 0);
    var has_dad = ("dad" in s && s.dad && s.dad.alts >= 0);

    if (has_mom || has_dad) {
      var mom_ok = has_mom && s.mom.affected && s.mom.het;
      var dad_ok = has_dad && s.dad.affected && s.dad.alts >= 1;
      // If at least one parent is present, one must be an affected carrier
      if (!mom_ok && !dad_ok) return false;
    }
    return true;
  }

  // AFFECTED MALES
  if (s.sex == "male" || s.sex == 1 || s.sex == "1") {
    // Affected males must carry the ALT allele
    if (s.alts < 1) return false;

    // Male-to-male transmission is impossible on X
    if ("dad" in s && s.dad && s.dad.alts >= 0) {
      if (s.dad.affected && s.dad.alts >= 1) return false;
    }

    // If mother is present & sequenced, she must be an affected carrier
    if ("mom" in s && s.mom && s.mom.alts >= 0) {
      if (!s.mom.affected && s.mom.alts > 0) return false;
    }

    return true;
  }

  // Fallback for unknown sex: must carry ALT
  return s.alts >= 1;
}


function hom_ref_parent(s) {
  return ("dad" in s) && s.dad.hom_ref && ("mom" in s) && s.mom.hom_ref;
}


function segregating_dominant(s) {
  // Delegate X-chromosome variants
  if (variant.CHROM === "chrX" || variant.CHROM === "X") { 
    return segregating_dominant_x(s); 
  }
  // Affected samples MUST be Heterozygous (0/1)
  if (s.affected) {
    return s.het;
  }
  // Unaffected samples MUST be Homozygous Reference (0/0)
  return s.hom_ref;
}

function segregating_recessive_x(s) {
  if (!s || s.unknown) return false;

  // 1. FEMALE LOGIC
  if (s.sex == "female" || s.sex == 2 || s.sex == "2") {
    if (s.affected) {
      // Affected female MUST be homozygous alternate (X^a X^a)
      if (!s.hom_alt && s.alts < 2) return false;

      // Transmission Check: Both parents (if present & sequenced)
      // - Mother must be at least a carrier (alts >= 1)
      if ("mom" in s && s.mom && s.mom.alts >= 0) {
        if (s.mom.alts < 1) return false;
      }
      // - Father MUST be affected and carry ALT (X^a Y)
      if ("dad" in s && s.dad && s.dad.alts >= 0) {
        if (!s.dad.affected || s.dad.alts < 1) return false;
      }

      return true;
    } else {
      // Unaffected female can be hom_ref (0/0), carrier het (0/1), or missing (./.)
      return s.het || s.hom_ref || s.alts < 0;
    }
  }

  // MALE LOGIC
  if (s.sex == "male" || s.sex == 1 || s.sex == "1") {
    if (s.affected) {
      // Affected male must carry the ALT allele (alts >= 1)
      if (s.alts < 1) return false;

      // Transmission Check: Mother (if present & sequenced)
      // Affected male inherits X exclusively from mother -> Mother cannot be 0/0
      if ("mom" in s && s.mom && s.mom.alts >= 0) {
        if (s.mom.alts < 1) return false;
      }

      return true;
    } else {
      // Unaffected male MUST NOT carry ALT allele (0/0 or missing ./.)
      return s.hom_ref || s.alts < 0;
    }
  }

  return false;
}


function segregating_recessive(s) {
  if (variant.CHROM == "chrX" || variant.CHROM == "X") {
    return segregating_recessive_x(s);
  }
  // AFFECTED INDIVIDUALS
  if (s.affected) {
    // Proband / affected sample MUST be homozygous alternate (1/1)
    if (!s.hom_alt) return false;

    // Transmission Check — Mother
    if ("mom" in s && s.mom && s.mom.alts >= 0) {
      if (s.mom.affected) {
        // Affected mother MUST be hom_alt (1/1)
        if (!s.mom.hom_alt) return false;
      } else {
        // Unaffected mother MUST be a carrier (0/1) — cannot be 0/0 or 1/1
        if (!s.mom.het) return false;
      }
    }
    // Transmission Check — Father
    if ("dad" in s && s.dad && s.dad.alts >= 0) {
      if (s.dad.affected) {
        // Affected father MUST be hom_alt (1/1)
        if (!s.dad.hom_alt) return false;
      } else {
        // Unaffected father MUST be a carrier (0/1) — cannot be 0/0 or 1/1
        if (!s.dad.het) return false;
      }
    }

    return true;
  }
  // UNAFFECTED INDIVIDUALS (Controls, siblings, non-proband relatives)
  // Unaffected members can be carrier (0/1), wild-type (0/0), or missing (./.)
  return s.het || s.hom_ref || s.alts < 0;
}


function parents_x_dn_or_homref(s) {
  // Mother must not carry any alternate alleles (alts == 0 or missing/alts == -1)
  if ("mom" in s && s.mom && s.mom.alts > 0) {
    return false;
  }
  // Father must not carry any alternate alleles
  if ("dad" in s && s.dad && s.dad.alts > 0) {
    return false;
  }
  return true;
}


function affected_het_leaf(s) {
  if ("mom"in s && !affected_het_leaf(s.mom)) { return false; }
  if ("dad" in s && !affected_het_leaf(s.dad)) { return false; }
  return true;
}


function segregating_denovo(s) {
  // Check unaffected samples
  if (!s.affected) { return s.hom_ref; }
  // Reject homozygous alt calls
  if (s.hom_alt) { return false; }
  // Check mother if present
  if ("mom" in s) {
    if (s.mom.affected || s.mom.alts > 0) { return false; }
  }
  // Check father if present
  if ("dad" in s) {
    if (s.dad.affected || s.dad.alts > 0) { return false; }
  }
  // Must be a complete trio to confirm a de novo call
  if (!("mom" in s) || !("dad" in s)) { return false; }
  // Affected proband must be heterozygous (0/1)
  return s.het;
}

function segregating_x_denovo(s) {
  if (variant.CHROM !== "chrX" && variant.CHROM !== "X") return false;
   
  // Both parents MUST be present and sequenced
  var has_mom = ("mom" in s && s.mom && s.mom.alts >= 0);
  var has_dad = ("dad" in s && s.dad && s.dad.alts >= 0);
  if (!has_mom || !has_dad) return false;

  // Both parents MUST be unaffected and wild-type (0/0)
  if (s.mom.affected || s.dad.affected) return false;
  if (s.mom.alts > 0 || s.dad.alts > 0) return false;

  // Proband check based on sex
  if (!s.affected) return false;

  if (s.sex == "female" || s.sex == 2 || s.sex == "2") {
    // Affected female de novo on X requires a heterozygous call (0/1)
    return s.het || s.alts == 1;
  }

  if (s.sex == "male" || s.sex == 1 || s.sex == "1") {
    // Affected male de novo on X requires hemizygous/ALT allele (>=1)
    return s.alts >= 1;
  }

  return false;
}
