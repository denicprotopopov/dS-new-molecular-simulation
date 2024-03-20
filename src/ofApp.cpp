#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){
    ofSetWindowTitle("particleMD");
    gui.setDefaultWidth(500);
    gui.setup();
    gui.add(targetTemperature.setup("equilibrium temperature", 25000.0, 1000.0, 1000000.0));
    gui.add(coupling.setup("coupling", 0.5, 0.1, 1.0));
    gui.add(applyThermostat.setup("Apply Thermostat", true));
    gui.add(particleAmount.setup("particle amount", 100, 2, 1000));
    gui.add(applyAmount.setup("Apply amount (reloads simulation for now)", false));
    gui.add(label.setup("Controls", "Press F to enter or exit fullscreen"));
    ofBackground(0);
    ofSetFrameRate(60);
    initializeParticles();
}

//--------------------------------------------------------------
void ofApp::update(){
    for (auto &particle : particles) {
            updateParticle(particle, timeStep);
        }

        resolveParticleCollisions();
    
        if(applyAmount){
            initializeParticles();
        }
    
        if(applyThermostat) {
            applyBerendsenThermostat();
        }
        for (auto &particle : particles) {
            checkWallCollisions(particle);
        }
}

//--------------------------------------------------------------
void ofApp::draw(){
    for (const auto &particle : particles) {
           // particle
           ofDrawCircle(particle.position, particle.radius);

            // trail
           ofSetLineWidth(8);
           float alphaStep = 255.0f / particle.trail.size();
           for (size_t i = 1; i < particle.trail.size(); i++) {
               float alpha = 255 - (i * alphaStep);
               ofSetColor(255, 255, 255, alpha);
               const glm::vec2 &currentPosition = particle.trail[i];
               const glm::vec2 &previousPosition = particle.trail[i-1];
               ofDrawLine(previousPosition, currentPosition);
           }
           // Reset drawing settings
           ofSetLineWidth(1); // Reset line width to default
           ofSetColor(255, 255, 255, 255); // Reset color
       }
    
    gui.draw();
}


glm::vec2 ofApp::computeForce(Particle &particle) {
    glm::vec2 mousePos(ofGetMouseX(), ofGetMouseY());
        return -0.2f * (particle.position - mousePos);
}


void ofApp::applyBerendsenThermostat() {
    float targetTemp = targetTemperature; // m_EquilibriumTemperature
    float tau = coupling; // m_BerendsenThermostatCoupling
    float kB = 8.314; // Boltzmann constant

    // temp calculation
    float currentTemperature = 0.0;
    float kineticEnergy = 0.0;
    for (const auto& particle : particles) {
        kineticEnergy += 0.5f * particle.mass * glm::length2(particle.velocity);
    }
    currentTemperature = kineticEnergy / (3.0 * kB);

    // velocity scaling factor
    //scaleFactor = (float)Math.Sqrt(m_EquilibriumTemperature / (m_BerendsenThermostatCoupling * temperature));
    float scaleFactor = sqrt(targetTemp / (tau * currentTemperature));
    for (auto& particle : particles) {
        particle.velocity *= scaleFactor;
    }
}

void ofApp::updateParticle(Particle &particle, float deltaTime) {
    glm::vec2 force = computeForce(particle);
    particle.velocity += deltaTime * force / particle.mass;
    particle.position += deltaTime * particle.velocity;

    particle.trail.push_front(particle.position);

        // removes old positions
        if (particle.trail.size() > particle.trailMaxLength) {
            particle.trail.pop_back();
        }
}

void ofApp::checkWallCollisions(Particle &particle) {
    if (particle.position.x < 0 || particle.position.x > ofGetWidth()) {
        particle.velocity.x *= -1.0f;
        particle.position.x = ofClamp(particle.position.x, 0, ofGetWidth());
    }
    if (particle.position.y < 0 || particle.position.y > ofGetHeight()) {
        particle.velocity.y *= -1.0f;
        particle.position.y = ofClamp(particle.position.y, 0, ofGetHeight());
    }
}

void ofApp::resolveParticleCollisions() {
    for (auto& particle : particles) {
            particle.hasCollided = false;
        }
    for (int i = 0; i < particles.size(); ++i) {
            for (int j = i + 1; j < particles.size(); ++j) {
                glm::vec2 delta = particles[j].position - particles[i].position;
                float dist = glm::length(delta);
                float totalRadius = particles[i].radius + particles[j].radius;

                // detecting collisions
                if (dist < totalRadius) {
                    float mi = particles[i].mass;
                    float mj = particles[j].mass;
                    glm::vec2 vi = particles[i].velocity;
                    glm::vec2 vj = particles[j].velocity;

                    //    calculate center of mass (COM) velocity
                    glm::vec2 vCOM = (mi * vi + mj * vj) / (mi + mj);
                    //    calculate velocity of particle i & j in COM frame
                    glm::vec2 viCOM = vi - vCOM;
                    glm::vec2 vjCOM = vj - vCOM;

                    // soft collisions
                    particles[i].velocity = vCOM - viCOM;
                    particles[j].velocity = vCOM - vjCOM;

                    // prefenting overlap
                    float overlap = 0.5f * (totalRadius - dist);
                    glm::vec2 correction = glm::normalize(delta) * overlap;
                    particles[i].position -= correction;
                    particles[j].position += correction;
                    particles[i].hasCollided = true;
                    particles[j].hasCollided = true;
                    
                }
            }
        }
}

void ofApp::initializeParticles() {
    particles.clear();

    for (int i = 0; i < particleAmount; i++) {
        Particle p;
        p.position = glm::vec2(ofRandomWidth(), ofRandomHeight());
        p.velocity = glm::vec2(ofRandom(-1.0f, 100.0f), ofRandom(-1.0f, 200.0f));
        p.radius = 5.0f;
        p.mass = 5.0f / 0.5f;  // radius/0.5f
        p.hasCollided = false;
        particles.push_back(p);
    }
}
//--------------------------------------------------------------
void ofApp::exit(){

}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
    if (key == 'f'){
        ofToggleFullscreen();
    }
}
