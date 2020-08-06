from flask import Flask
import  GenomeAnalysis
app = Flask(__name__)


@app.route("/fhirgenome/genomesequence/<genomeId>")
def getGenome(genomeId):
    return GenomeAnalysis.returnGenome(genomeId)

@app.route("/fhirgenome/oric/<genomeId>")
def getOric(genomeId):
    return GenomeAnalysis.returnOric(genomeId)

@app.route("/fhirgenome/transcription/<genomeId>")
def transcription(genomeId):
    return GenomeAnalysis.transcription(genomeId)

@app.route("/fhirgenome/translation/<genomeId>")
def translation(genomeId):
    return GenomeAnalysis.translation(genomeId)

@app.route("/fhirgenome/polyPeptides/<genomeId>")
def aminoAcids(genomeId):
    return GenomeAnalysis.amino_acids(genomeId).to_dict()

@app.route("/fhirgenome/largestPolyPeptides/<genomeId>")
def largestAminoAcids(genomeId):
    return GenomeAnalysis.top_ten_largest_amino_acids(genomeId).to_dict()

@app.route("/fhirgenome/topAminoAcids/<genomeId>")
def topAminoAcids(genomeId):
    return GenomeAnalysis.top_amino_acids(genomeId)

if __name__ == '__main__':
    app.run(debug=True,host="0.0.0.0",port="5001")