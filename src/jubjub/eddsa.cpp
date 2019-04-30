// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "jubjub/eddsa.hpp"
#include "utils.hpp"

namespace ethsnarks {

namespace jubjub {

//Private key
//An EdDSA private key is a b-bit string k which should be chosen uniformly at random.

//Public key
//The corresponding public key is a curve point AE(Fq), encoded in bbits. A=sB, where s=H0,.........,Hb-1(k) is the least significant bbits of H(k) interpreted as integer in little-endian
    
//Signature
//An EdDSA signature on a message m by public key A is the pair (R,S), encoded in 2bbits, of a curve point RE(Fq) and an integer 0<S<l satisfying the verification equation, where 
//R=rB for r=H(Hb,......,H2b-1(k),M) and S=r+H(R,A,M)s mod l. 

//Verification:
//2cSB=2c(r+H(R,A,M)s)B= 2crB+2cH(R,A,M)sB= 2cR+2cH(R,A,M)A

//we define the protoboard and the variables we need in the hash function H(R,A,M) which will be needed afterwards for the signature S
EdDSA_HashRAM_gadget::EdDSA_HashRAM_gadget(
    ProtoboardT& in_pb,
    const Params& in_params,
    const VariablePointT& in_R,
    const VariablePointT& in_A,
    const VariableArrayT& in_M,
    const std::string& annotation_prefix
) :
    GadgetT(in_pb, annotation_prefix),

    // Convert X & Y coords to bits for hash function
    m_R_x_bits(in_pb, in_R.x, FMT(this->annotation_prefix, ".R_x_bits")),
    m_A_x_bits(in_pb, in_A.x, FMT(this->annotation_prefix, ".A_x_bits")),

    // Prefix the message with R and A.
    m_RAM_bits(flatten({
        m_R_x_bits.result(),
        m_A_x_bits.result(),
        in_M,
    })),

    m_hash_RAM(in_pb, in_params, "EdDSA_Verify.RAM", m_RAM_bits, FMT(this->annotation_prefix, ".hash_RAM"))
{
}

//generate_r1cs_constraints() function adds the R1CS constraints corresponding to the circuits defined earlier.
void EdDSA_HashRAM_gadget::generate_r1cs_constraints()
{
    m_R_x_bits.generate_r1cs_constraints();
    m_A_x_bits.generate_r1cs_constraints();
    m_hash_RAM.generate_r1cs_constraints();
}

//generate_r1cs_witness() function assumes that we've already set the public value and the witness value. It then computes the inferred witness values for the intermediate variables 


void EdDSA_HashRAM_gadget::generate_r1cs_witness()
{
    m_R_x_bits.generate_r1cs_witness();
    m_A_x_bits.generate_r1cs_witness();
    m_hash_RAM.generate_r1cs_witness();
}


const VariableArrayT& EdDSA_HashRAM_gadget::result()
{
    return m_hash_RAM.result();
}


// --------------------------------------------------------------------
//  One of the parameters of the EdDSA algorithm is the "prehash" function.  This may be the identity function, resulting in an
   //algorithm called PureEdDSA, or a collision-resistant hash function such as SHA-512, resulting in an algorithm called HashEdDSA.

//we define the protoboard and the variables we need in the signature S=r+H(R,A,M)s mod l for the pureEddsa algorithm
    // In PureEdDSA algorithm, there is no message compression function (M = H'(m)) to compress m.

// Ed25519 and Ed448 are PureEdDSA variants 


PureEdDSA::PureEdDSA(
    ProtoboardT& in_pb,
    const Params& in_params,
    const EdwardsPoint& in_base,    // B
    const VariablePointT& in_A,     // A=s.B
    const VariablePointT& in_R,     // R=rB for r=H(Hb,......,H2b-1(k),m) 
    const VariableArrayT& in_s,     // s=H0,.........,Hb-1(k) is the least significant bbits of H(k) interpreted as integer in little-endian
    const VariableArrayT& in_msg,   // m
    const std::string& annotation_prefix
) :
    GadgetT(in_pb, annotation_prefix),

    // IsValid(R) to verify that R is a valid point on the curve
  

    m_validator_R(in_pb, in_params, in_R.x, in_R.y, FMT(this->annotation_prefix, ".validator_R")),

    // lhs = ScalarMult(B, s) where A=sB and s=H0,.........,Hb-1(k) is the least significant bbits of H(k) interpreted as integer in little-endian


    m_lhs(in_pb, in_params, in_base.x, in_base.y, in_s, FMT(this->annotation_prefix, ".lhs")),

    // hash_RAM = H(R, A, M)
    m_hash_RAM(in_pb, in_params, in_R, in_A, in_msg, FMT(this->annotation_prefix, ".hash_RAM")),

    // At = ScalarMult(A,hash_RAM)
    // since Base point B has order l so B mod l = 1 which means A mod l= s.B mod l = s
    // Which means H(R,A,M)s = H(R,A,M)A

    m_At(in_pb, in_params, in_A.x, in_A.y, m_hash_RAM.result(), FMT(this->annotation_prefix, ".At = A * hash_RAM")),

    // rhs = PointAdd(R, At) which represents the signature
    
    m_rhs(in_pb, in_params, in_R.x, in_R.y, m_At.result_x(), m_At.result_y(), FMT(this->annotation_prefix, ".rhs"))
{ }


void PureEdDSA::generate_r1cs_constraints()
{
    m_validator_R.generate_r1cs_constraints();
    m_lhs.generate_r1cs_constraints();
    m_hash_RAM.generate_r1cs_constraints();
    m_At.generate_r1cs_constraints();
    m_rhs.generate_r1cs_constraints();

    // Verify the two points are equal
    this->pb.add_r1cs_constraint(
        ConstraintT(m_lhs.result_x(), FieldT::one(), m_rhs.result_x()),
        FMT(this->annotation_prefix, " lhs.x == rhs.x"));

    this->pb.add_r1cs_constraint(
        ConstraintT(m_lhs.result_y(), FieldT::one(), m_rhs.result_y()),
        FMT(this->annotation_prefix, " lhs.y == rhs.y"));
}


void PureEdDSA::generate_r1cs_witness()
{
    m_validator_R.generate_r1cs_witness();
    m_lhs.generate_r1cs_witness();
    m_hash_RAM.generate_r1cs_witness();
    m_At.generate_r1cs_witness();
    m_rhs.generate_r1cs_witness();
}


// --------------------------------------------------------------------

//HashEdDSA or edDSA requires an extra step to compress the message before hashing, this is: M = H'(m)
EdDSA::EdDSA(
    ProtoboardT& in_pb,
    const Params& in_params,
    const EdwardsPoint& in_base,    // B
    const VariablePointT& in_A,     // A
    const VariablePointT& in_R,     // R
    const VariableArrayT& in_s,     // s
    const VariableArrayT& in_msg,   // m
    const std::string& annotation_prefix
) :
    // M = H(m)
    m_msg_hashed(in_pb, in_params, "EdDSA_Verify.M", in_msg, FMT(annotation_prefix, ".msg_hashed")),

    m_verifier(in_pb, in_params, in_base, in_A, in_R, in_s, m_msg_hashed.result(), annotation_prefix)
{ }


void EdDSA::generate_r1cs_constraints()
{
    m_msg_hashed.generate_r1cs_constraints();
    m_verifier.generate_r1cs_constraints();
}


void EdDSA::generate_r1cs_witness()
{
    m_msg_hashed.generate_r1cs_witness();
    m_verifier.generate_r1cs_witness();
}


// namespace jubjub
}

// namespace ethsnarks
}
